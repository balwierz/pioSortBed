#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cstdarg>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <execution>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <time.h>
#include <unordered_map>
#include <unordered_set>
#include <errno.h>
#include <tbb/global_control.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "CLI11.hpp"
#ifdef WITH_BAM
#include <htslib/sam.h>
#include <htslib/hts.h>
#endif
#ifdef WITH_RANS
#include <htscodecs/rANS_static4x16.h>
#endif
#include <lz4.h>
#include <zstd.h>

// Fallback for manual compilation without -DVERSION_STRING; the Makefile
// always injects the real version, so this only matters for one-off builds.
#ifndef VERSION_STRING
#define VERSION_STRING "manual-build"
#endif

/* Copyright: Piotr Balwierz */

using namespace std;

// ============================================================================
// CORE TYPES & ALLOCATOR
// ----------------------------------------------------------------------------
// seqread:    one parsed input line (24 B; lives in reads[]).
// chrInfoT:   per-chromosome bookkeeping (linked-list head).
// ChrNameMap: small flat-vector map keyed by chromosome name.
// Arena:      bump allocator for line tails / weight strings on stdin/gzip.
// Constants:  weight-field cap.
// ============================================================================

// Field layout minimizes padding: pointer first, then ints, then char.
// 24 bytes per read on 64-bit (8 + 4 + 4 + 4 + 1 + 3 padding).
//
// The `next` field (linked lists during parsing) and `chrIdx` field (classic
// sort path integer chromosome index) share the same 4 bytes via a union —
// the two sort paths never need both simultaneously. Both are unsigned 32-bit:
// `next` is an index into reads[], 0 = end-of-list sentinel; `chrIdx` is a
// sequential chromosome index assigned after parsing.
class seqread
{
	public:
	char* line; // keeps the full line; in case of collapse option just the weight as string.
	int beg;
	int end;
	union { uint32_t next; uint32_t chrIdx; };
	char str;
};

class chrInfoT
{
	public:
	uint32_t lastRead;	// head of the linked list of read indices for this chromosome
						// (0 = no reads yet / end-of-list sentinel)
	uint32_t idx;		// sequential chromosome index, assigned after parsing & alphabetical sort
	chrInfoT() : lastRead(0), idx(0) {}
};

// Small flat-vector "map" keyed by chromosome name.
//
// Real BED files have <100 chromosomes, so the linear scan beats the
// std::unordered_map<string, V> alternative by avoiding both the per-lookup
// hash + bucket probe and (more importantly) the std::string construction
// the unordered_map forces on every miss. The reserved capacity (1024) is
// far above the practical chromosome count; we never resize it, so cached
// iterators stay valid across emplaces.
template<typename Value>
class ChrNameMap
{
public:
	using Entry = std::pair<std::string, Value>;
	using iterator = typename std::vector<Entry>::iterator;
	using const_iterator = typename std::vector<Entry>::const_iterator;

	std::vector<Entry> entries;

	ChrNameMap() { entries.reserve(1024); }

	iterator       begin()       { return entries.begin(); }
	iterator       end()         { return entries.end();   }
	const_iterator begin() const { return entries.begin(); }
	const_iterator end()   const { return entries.end();   }
	size_t         size()  const { return entries.size();  }
	void           clear()       { entries.clear();        }

	// Find by (ptr, len) — avoids constructing a std::string at the call site.
	iterator findByPtr(const char* name, int len)
	{
		for(auto it = entries.begin(); it != entries.end(); ++it)
			if((int)it->first.size() == len && memcmp(it->first.data(), name, len) == 0)
				return it;
		return entries.end();
	}
	iterator find(const std::string& name)
	{
		return findByPtr(name.data(), (int)name.size());
	}

	// Emplace a new entry only if (ptr, len) isn't already present.
	// Returns (iterator, true_if_inserted), matching std::unordered_map::try_emplace.
	std::pair<iterator, bool> tryEmplaceByPtr(const char* name, int len)
	{
		iterator it = findByPtr(name, len);
		if(it != entries.end()) return {it, false};
		entries.emplace_back(std::string(name, len), Value{});
		return {entries.end() - 1, true};
	}
	std::pair<iterator, bool> try_emplace(const std::string& name)
	{
		return tryEmplaceByPtr(name.data(), (int)name.size());
	}
	std::pair<iterator, bool> insert(const std::pair<std::string, Value>& kv)
	{
		iterator it = find(kv.first);
		if(it != entries.end()) return {it, false};
		entries.push_back(kv);
		return {entries.end() - 1, true};
	}

	size_t erase(const std::string& name)
	{
		iterator it = find(name);
		if(it == entries.end()) return 0;
		entries.erase(it);
		return 1;
	}
};

typedef ChrNameMap<chrInfoT> string2chrInfoT;

// Stack buffer used for the BED weight field (col 4) on the --collapse path.
// Real weights are short numeric strings; 256 B is generous. On overflow,
// copyField returns 0 and the parser errors with a clear message rather than
// silently truncating.
const int kWeightBufSize = 256;


// Arena allocator for line/weight strings.
// Used for stdin/gzip input (all lines) and for collapse mode (weight strings).
// Not used for file input in non-collapse mode (lines point into mmap).
struct ArenaChunk
{
	char*       data;
	size_t      used;
	size_t      capacity;
	ArenaChunk* next;

	ArenaChunk(size_t cap) : used(0), capacity(cap), next(NULL)
	{
		data = (char*) malloc(cap);
		if (!data) { fprintf(stderr, "Arena: out of memory\n"); exit(1); }
	}
	~ArenaChunk() { free(data); }
};

struct Arena
{
	ArenaChunk* head;
	size_t      chunkSize;

	Arena(size_t initialCap, size_t chunkSz) : chunkSize(chunkSz)
	{
		head = new ArenaChunk(initialCap);
	}
	~Arena()
	{
		ArenaChunk* c = head;
		while (c) { ArenaChunk* nx = c->next; delete c; c = nx; }
	}

	char* alloc(const char* src, size_t len)
	{
		ArenaChunk* c = head;
		while (c->used + len > c->capacity)
		{
			if (c->next == NULL)
			{
				size_t newCap = (len > chunkSize) ? len : chunkSize;
				c->next = new ArenaChunk(newCap);
			}
			c = c->next;
		}
		char* ptr = c->data + c->used;
		memcpy(ptr, src, len);
		c->used += len;
		return ptr;
	}
};


// ============================================================================
// PARSERS — FIELDS, BED LINES, COLLAPSE-MODE WEIGHTS
// ----------------------------------------------------------------------------
// Field-level (skipField/parseUInt/copyField) and line-level
// (parseBedLine3/parseBedLineFull) primitives, all `inline` for the hot path.
// Plus parseWeight/sumWeights* which the bucket-sort and classic-sort emit
// paths use to consume `--collapse` weight strings. Hand-written; faster
// than sscanf.
// ============================================================================

// Hand-written field parsers (faster than sscanf).
static inline bool skipField(const char** pp)
{
	const char* p = *pp;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	if (*p == '\t' || *p == ' ') { *pp = p + 1; return true; }
	*pp = p;
	return false;
}

static inline int parseUInt(const char** pp)
{
	const char* p = *pp;
	if (*p < '0' || *p > '9') return -1;
	int val = 0;
	while (*p >= '0' && *p <= '9') { val = val * 10 + (*p - '0'); p++; }
	*pp = p;
	return val;
}

// Copy a field into dst (NUL-terminated). Returns the number of bytes
// written, or -1 if the field would not fit in dstMax (caller errors).
// Always advances *pp past the trailing separator (or to EOL/EOF).
static inline int copyField(const char** pp, char* dst, int dstMax)
{
	const char* p = *pp;
	int n = 0;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n')
	{
		if (n >= dstMax - 1) {
			// Overflow: skip the rest of the field so *pp lands on the
			// separator, but signal failure to the caller.
			while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
			if (*p == '\t' || *p == ' ') p++;
			*pp = p;
			return -1;
		}
		dst[n++] = *p;
		p++;
	}
	dst[n] = '\0';
	if (*p == '\t' || *p == ' ') p++;
	*pp = p;
	return n;
}

// Minimal BED parser: chr (as pointer+length), beg, end, tailPtr (fields 4+).
static inline int parseBedLine3(const char* buf,
                                const char** chrPtr, int* chrLen,
                                int* beg, int* end,
                                const char** tailPtr)
{
	const char* p = buf;
	*chrPtr = p;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	*chrLen = (int)(p - *chrPtr);
	if (*chrLen == 0) return 0;
	if (*p == '\t' || *p == ' ') p++; else return 0;

	int v = parseUInt(&p);
	if (v < 0) return 1;
	*beg = v;
	if (*p == '\t' || *p == ' ') p++;

	v = parseUInt(&p);
	if (v < 0) return 2;
	*end = v;
	*tailPtr = p;
	return 3;
}

// Full BED parser: chr (col 1), beg (col 2), end (col 3), weight (col 5,
// the score), strand (col 6). Col 4 (name) is skipped — never stored.
static inline int parseBedLineFull(const char* buf,
                                   const char** chrPtr, int* chrLen,
                                   int* beg, int* end,
                                   char* weightBuf, int weightMax,
                                   char* strandChar,
                                   const char** tailPtr)
{
	const char* p = buf;
	*chrPtr = p;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	*chrLen = (int)(p - *chrPtr);
	if (*chrLen == 0) return 0;
	if (*p == '\t' || *p == ' ') p++; else return 0;

	int v = parseUInt(&p);
	if (v < 0) return 1;
	*beg = v;
	if (*p == '\t' || *p == ' ') p++;

	v = parseUInt(&p);
	if (v < 0) return 2;
	*end = v;
	*tailPtr = p;

	if (*p != '\t' && *p != ' ') return 3;
	p++;

	// col 4 (name) — skip it
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	if (*p != '\t' && *p != ' ') return 3;
	p++;

	// col 5 (score, used as weight by --collapse): stack buffer, over-long
	// values rejected by copyField (BED scores are spec'd 0..1000 = at most
	// 4 chars; 255 B is a generous ceiling for non-spec extensions).
	int n = copyField(&p, weightBuf, weightMax);
	if (n == 0) return 3;
	if (n < 0) {
		fprintf(stderr, "Error: BED score field (column 5) exceeds %d bytes in: %s\n",
		        weightMax - 1, buf);
		exit(1);
	}

	// col 6 (strand) — copyField already consumed the tab before this field
	if (!*p || *p == '\r' || *p == '\n') return 4;
	*strandChar = *p;
	return 5;
}

// Parse a weight string written into reads[idx].line (collapse mode) using
// strtof (faster + locale-independent compared to sscanf). Aborts on malformed input.
static inline float parseWeight(const seqread* reads, uint32_t idx)
{
	char* endptr;
	float w = strtof(reads[idx].line, &endptr);
	if(endptr == reads[idx].line)
	{
		cerr << "Malformed weight " << reads[idx].line << endl;
		exit(1);
	}
	return w;
}

// ============================================================================
// LOW-MEMORY SSD PATH (--low-mem-ssd)
// ----------------------------------------------------------------------------
// Two-pass file-only sort. Pass 1: parse mmap once, build a 16 B-per-line
// node table with per-chromosome linked lists threaded through it. Pass 2:
// for each chromosome, walk its list, sort the small per-chrom vector, emit
// the lines via fwrite directly from mmap pointers. Pass 2 runs chromosomes
// in parallel; output is serialised through a producer-consumer barrier.
//
// Owns: lowMemNode (the per-line index entry), naturalChrLess (the
// chromosome ordering helper), ChromBuf (the per-chrom output buffer),
// lowMemSortMmap (the whole pipeline).
// ============================================================================

// One node per line in --low-mem-ssd's pass-1 index.
//
//   off  — byte offset of the line's start within the mmap (size_t so files
//          can exceed 4 GB; ~100M+ reads).
//   next — index of the next node in this chromosome's linked list (0 = end).
//   beg  — pre-parsed BED start coordinate. The default --sort s mode and
//          --collapse only need beg, so caching it lets those modes sort
//          and emit without re-parsing the line. --sort b re-parses end
//          on demand (cheap: one tab-scan + atoi per read in pass 2).
//          --sort 5 already re-parses for strand, and pulls end out of
//          that same parse for free.
//
// 16 bytes / node — exactly aligned, no padding. 50M reads cost ~0.8 GB of
// node table; 200M reads cost ~3.2 GB (down from ~4.8 GB before end was
// dropped).
struct lowMemNode
{
	size_t off;
	uint32_t next;
	int beg;
};

// Natural (version-sort) comparison for chromosome names.
// Embedded digit runs are compared numerically: "chr2" < "chr10".
static bool naturalChrLess(const std::string& a, const std::string& b)
{
	const char* pa = a.c_str();
	const char* pb = b.c_str();
	while (*pa && *pb)
	{
		if (isdigit((unsigned char)*pa) && isdigit((unsigned char)*pb))
		{
			// skip leading zeros so "01" == "1" numerically
			while (*pa == '0' && isdigit((unsigned char)*(pa+1))) ++pa;
			while (*pb == '0' && isdigit((unsigned char)*(pb+1))) ++pb;
			// measure digit-run lengths — longer run = larger number
			const char* ea = pa; while (isdigit((unsigned char)*ea)) ++ea;
			const char* eb = pb; while (isdigit((unsigned char)*eb)) ++eb;
			int la = (int)(ea - pa), lb = (int)(eb - pb);
			if (la != lb) return la < lb;
			int c = strncmp(pa, pb, la);
			if (c) return c < 0;
			pa = ea; pb = eb;
		}
		else
		{
			if (*pa != *pb) return (unsigned char)*pa < (unsigned char)*pb;
			++pa; ++pb;
		}
	}
	return (unsigned char)*pa < (unsigned char)*pb;
}

// Per-chromosome output buffer used by --low-mem-ssd's pass 2.
//
// Two modes:
//   flushTo == NULL  — accumulate everything into `buf`. Used by the multi-
//                      thread path so the producer-consumer barrier can flush
//                      whole chromosomes to stdout in alphabetical order.
//                      Pre-size the buffer (see avgLineBytes estimate at the
//                      call site) to skip the early doublings.
//   flushTo != NULL  — small reusable scratch buffer; when it fills we
//                      fwrite_unlocked to flushTo and reset pos. Used by the
//                      single-thread path to bound peak RAM at ~64 KB instead
//                      of ~per-chromosome.
//
// Replaces the previous `open_memstream` approach: skips the stdio FILE state
// machine and (in the buffered case) avoids the realloc churn that the
// memstream's exponential doubling caused on multi-GB chromosomes.
struct ChromBuf
{
	char* buf = nullptr;
	size_t pos = 0;
	size_t cap = 0;
	FILE* flushTo = nullptr;

	void grow_or_flush(size_t need)
	{
		if(flushTo)
		{
			if(pos > 0)
			{
				fwrite_unlocked(buf, 1, pos, flushTo);
				pos = 0;
			}
			if(need <= cap) return;
		}
		size_t newCap = cap ? cap : 65536;
		while(newCap < pos + need) newCap *= 2;
		char* nb = (char*) realloc(buf, newCap);
		if(!nb) { perror("realloc ChromBuf"); exit(1); }
		buf = nb;
		cap = newCap;
	}

	inline void writeLineFromMmap(const char* line)
	{
		size_t len = strlen(line);
		if(pos + len + 1 > cap) grow_or_flush(len + 1);
		memcpy(buf + pos, line, len);
		buf[pos + len] = '\n';
		pos += len + 1;
	}

	void writeFmt(const char* fmt, ...)
	{
		char tmp[512];
		va_list ap;
		va_start(ap, fmt);
		int n = vsnprintf(tmp, sizeof(tmp), fmt, ap);
		va_end(ap);
		if(n < 0) return;
		if((size_t)n < sizeof(tmp))
		{
			if(pos + (size_t)n > cap) grow_or_flush((size_t)n);
			memcpy(buf + pos, tmp, (size_t)n);
			pos += (size_t)n;
			return;
		}
		// Output exceeded 512 chars (rare). Format a second time directly into buf.
		if(pos + (size_t)n + 1 > cap) grow_or_flush((size_t)n + 1);
		va_start(ap, fmt);
		int n2 = vsnprintf(buf + pos, cap - pos, fmt, ap);
		va_end(ap);
		if(n2 < 0) return;
		pos += (size_t)n2;
	}

	void finalize()
	{
		if(flushTo && pos > 0)
		{
			fwrite_unlocked(buf, 1, pos, flushTo);
			pos = 0;
		}
	}
};

// Low-memory file mode (SSD-friendly):
// Pass 1: walk mmap once, build per-chromosome linked lists of line offsets.
// Pass 2: process one chromosome at a time (parse, sort, print), so peak RAM
// depends on the largest chromosome chunk, not the whole file.
static int lowMemSortMmap(char* mmapBase, size_t mmapSize,
                          int fCollapse, char sortMode, int numThreads,
                          bool naturalSort, size_t maxMemBytes, bool verbose)
{
	time_t tstart, tend;
	if(verbose) time(&tstart);

	string2chrInfoT chrInfo;
	lowMemNode* nodes = NULL;
	size_t nodeCount = 1;  // index 0 reserved as end-of-list sentinel
	// Hard cap: node indices are stored in lowMemNode::next (uint32_t), so the
	// usable index range is 1..UINT32_MAX-1 → at most ~4.29 B reads.

	// Pre-pass: emit any leading BED header lines (#, track, browser) to stdout
	// and advance bodyStart past them. Headers are leading-only by convention
	// (matches the regular --threads path); mid-file headers would now be
	// treated as data and trip the parser.
	size_t bodyStart = 0;
	while(bodyStart < mmapSize)
	{
		char* lineStart = mmapBase + bodyStart;
		char* nl = (char*) memchr(lineStart, '\n', mmapSize - bodyStart);
		size_t lineLen = (size_t)((nl ? nl : mmapBase + mmapSize) - lineStart);
		bool isHeader = (lineLen > 0 && lineStart[0] == '#') ||
		                (lineLen >= 6 && memcmp(lineStart, "track ",   6) == 0) ||
		                (lineLen >= 8 && memcmp(lineStart, "browser ", 8) == 0);
		if(!isHeader) break;
		fwrite_unlocked(lineStart, 1, lineLen, stdout);
		fputc_unlocked('\n', stdout);
		bodyStart = (nl ? (size_t)(nl - mmapBase) + 1 : mmapSize);
	}

	// Decide whether to use the parallel chunk parser. Below ~256 KB body the
	// scheduling overhead exceeds the parsing savings, so stay serial.
	int N = numThreads;
	if(N < 1) N = 1;
	{
		size_t bytesPerChunk = 256 * 1024;
		size_t bodySize = mmapSize - bodyStart;
		int byBytes = (int)((bodySize + bytesPerChunk - 1) / bytesPerChunk);
		if(byBytes < 1) byBytes = 1;
		if(N > byBytes) N = byBytes;
	}

	if(N > 1)
	{
		// ---- Pass 1, parallel ----
		// Newline-aligned chunk boundaries within [bodyStart, mmapSize).
		std::vector<size_t> chunkStart(N + 1);
		chunkStart[0] = bodyStart;
		chunkStart[N] = mmapSize;
		for(int i = 1; i < N; i++)
		{
			size_t guess = bodyStart + ((mmapSize - bodyStart) / (size_t)N) * (size_t)i;
			while(guess < mmapSize && mmapBase[guess] != '\n') guess++;
			if(guess < mmapSize) guess++;
			chunkStart[i] = guess;
		}

		// Pass 1a: count newlines per chunk in parallel so we can size nodes[] exactly.
		std::vector<size_t> chunkCount(N, 0);
		std::vector<int> indices(N);
		for(int i = 0; i < N; i++) indices[i] = i;
		std::for_each(std::execution::par, indices.begin(), indices.end(),
			[&](int i) {
				size_t a = chunkStart[i], b = chunkStart[i + 1];
				size_t c = 0;
				const char* p = mmapBase + a;
				const char* pe = mmapBase + b;
				while(p < pe)
				{
					const char* q = (const char*) memchr(p, '\n', (size_t)(pe - p));
					if(!q) break;
					c++;
					p = q + 1;
				}
				chunkCount[i] = c;
			});
		if(mmapSize > bodyStart && mmapBase[mmapSize - 1] != '\n')
			chunkCount[N - 1]++;

		// Prefix-sum -> per-chunk base offset in nodes[]. Slot 0 is the sentinel.
		std::vector<size_t> chunkBase(N + 1);
		chunkBase[0] = 0;
		for(int i = 0; i < N; i++) chunkBase[i + 1] = chunkBase[i] + chunkCount[i];
		size_t totalNodes = chunkBase[N];

		// Indices into nodes[] are stored in lowMemNode::next (uint32_t), so
		// we can't address more than UINT32_MAX-1 reads. 1 is the sentinel slot.
		if(totalNodes >= (size_t)UINT32_MAX)
		{
			cerr << "Error: " << totalNodes << " reads exceeds the 4.29 B limit "
			     << "imposed by the 32-bit per-line index. Split your input "
			     << "or sort chromosome-by-chromosome." << endl;
			return 1;
		}

		nodes = (lowMemNode*) malloc((totalNodes + 1) * sizeof(lowMemNode));
		if(!nodes)
		{
			cerr << "Error: out of memory allocating " << totalNodes
			     << "-entry low-memory offset table" << endl;
			return 1;
		}

		// Per-chunk per-chrom partials, indexed by global node index.
		struct LowMemChunkPartial { uint32_t head; uint32_t tail; };
		using LowMemChunkMap = ChrNameMap<LowMemChunkPartial>;
		std::vector<LowMemChunkMap> chunkChr(N);

		// Pass 1b: parse chunks in parallel into their slice of nodes[].
		std::for_each(std::execution::par, indices.begin(), indices.end(),
			[&](int i) {
				char* end = mmapBase + chunkStart[i + 1];
				char* mmapCur = mmapBase + chunkStart[i];
				uint32_t idx = (uint32_t)(chunkBase[i] + 1);
				LowMemChunkMap& chrMap = chunkChr[i];
				auto thisIt = chrMap.end();

				while(mmapCur < end)
				{
					char* linePtr = mmapCur;
					char* nl = (char*) memchr(mmapCur, '\n', (size_t)(end - mmapCur));
					if(nl)
					{
						if(nl > mmapCur && *(nl - 1) == '\r') *(nl - 1) = '\0';
						*nl = '\0';
						mmapCur = nl + 1;
					}
					else
					{
						mmapCur = end;
					}

					int beg = 0, lineEnd = 0;
					const char* chrPtr;
					int chrLen;
					const char* tailPtr = "";
					int numArgsRead = parseBedLine3(linePtr, &chrPtr, &chrLen,
					                                &beg, &lineEnd, &tailPtr);
					if(numArgsRead < 3)
					{
						cerr << "Error in parsing line: " << linePtr << endl;
						exit(1);
					}
					if(beg < 0)
					{
						cerr << "Error: negative coordinates in the bed file\n"
						     << linePtr << endl;
						exit(1);
					}

					(void)lineEnd;  // parsed for validation only; end isn't cached
					nodes[idx].off = (size_t)(linePtr - mmapBase);
					nodes[idx].beg = beg;

					if(thisIt != chrMap.end() &&
					   (int)thisIt->first.size() == chrLen &&
					   memcmp(thisIt->first.data(), chrPtr, chrLen) == 0)
					{
						nodes[idx].next = thisIt->second.head;
						thisIt->second.head = idx;
					}
					else
					{
						auto ins = chrMap.tryEmplaceByPtr(chrPtr, chrLen);
						thisIt = ins.first;
						if(ins.second)
						{
							nodes[idx].next = 0;
							thisIt->second.head = idx;
							thisIt->second.tail = idx;
						}
						else
						{
							nodes[idx].next = thisIt->second.head;
							thisIt->second.head = idx;
						}
					}
					idx++;
				}
			});

		// Merge per-chunk partials into the global chrInfo. For each chunk t,
		// for each chromosome it touched, prepend its list to the global list:
		// patch chunk's tail.next -> previous global head, then head = chunk's head.
		for(int t = 0; t < N; t++)
		{
			for(auto& kv : chunkChr[t])
			{
				const std::string& name = kv.first;
				const LowMemChunkPartial& partial = kv.second;
				auto ins = chrInfo.tryEmplaceByPtr(name.data(), (int)name.size());
				chrInfoT& info = ins.first->second;
				if(ins.second || info.lastRead == 0)
					info.lastRead = partial.head;
				else
				{
					nodes[partial.tail].next = info.lastRead;
					info.lastRead = partial.head;
				}
			}
		}

		nodeCount = totalNodes + 1;
	}
	else
	{
		// ---- Pass 1, serial (existing realloc-grow path) ----
		auto seedIns = chrInfo.insert(make_pair(string(""), chrInfoT()));
		auto thisChrIt = seedIns.first;

		size_t nodeCap = 1024;
		nodes = (lowMemNode*) malloc(nodeCap * sizeof(lowMemNode));
		if(!nodes)
		{
			cerr << "Error: out of memory allocating low-memory offset table" << endl;
			return 1;
		}

		char* mmapCur = mmapBase + bodyStart;
		char* mmapLim = mmapBase + mmapSize;
		while(mmapCur < mmapLim)
		{
			char* linePtr = mmapCur;
			char* nl = (char*) memchr(mmapCur, '\n', (size_t)(mmapLim - mmapCur));
			if(nl)
			{
				if(nl > mmapCur && *(nl - 1) == '\r') *(nl - 1) = '\0';
				*nl = '\0';
				mmapCur = nl + 1;
			}
			else
			{
				mmapCur = mmapLim;
			}

			const char* chrPtr;
			int chrLen;
			int beg = 0;
			int end = 0;
			const char* tailPtr = "";
			int numArgsRead = parseBedLine3(linePtr, &chrPtr, &chrLen, &beg, &end, &tailPtr);

			if(numArgsRead < 3)
			{
				cerr << "Error in parsing line: " << linePtr << endl
					<< "Perhaps this line is malformed?" << endl;
				free(nodes);
				return 1;
			}
			if(beg < 0)
			{
				cerr << "Error: negative coordinates in the bed file\n" << linePtr << endl;
				free(nodes);
				return 1;
			}

			if((int)thisChrIt->first.size() != chrLen ||
			   memcmp(thisChrIt->first.data(), chrPtr, chrLen) != 0)
			{
				// tryEmplaceByPtr avoids constructing a std::string when the chrom
				// already exists; we only allocate one if we actually insert.
				auto ins = chrInfo.tryEmplaceByPtr(chrPtr, chrLen);
				thisChrIt = ins.first;
			}

			if(nodeCount == nodeCap)
			{
				nodeCap *= 2;
				lowMemNode* tmp = (lowMemNode*) realloc(nodes, nodeCap * sizeof(lowMemNode));
				if(!tmp)
				{
					cerr << "Error: out of memory growing low-memory offset table" << endl;
					free(nodes);
					return 1;
				}
				nodes = tmp;
			}

			// nodeCount is the index we're about to write to nodes[]; it's also
			// the value we'll store as a uint32_t into chrInfoT::lastRead, so
			// it must fit in uint32_t. UINT32_MAX is reserved as well-out-of-range;
			// 0 is the end-of-list sentinel; usable range is 1..UINT32_MAX-1.
			if(nodeCount >= (size_t)UINT32_MAX)
			{
				cerr << "Error: " << nodeCount << " reads exceeds the 4.29 B limit "
				     << "imposed by the 32-bit per-line index." << endl;
				free(nodes);
				return 1;
			}

			(void)end;  // parsed for validation only; end isn't cached in lowMemNode
			nodes[nodeCount].off = (size_t)(linePtr - mmapBase);
			nodes[nodeCount].next = thisChrIt->second.lastRead;
			nodes[nodeCount].beg = beg;
			thisChrIt->second.lastRead = (uint32_t)nodeCount;
			nodeCount++;
		}

		chrInfo.erase("");
	}
	vector<string> chroms;
	chroms.reserve(chrInfo.size());
	for(string2chrInfoT::iterator it = chrInfo.begin(); it != chrInfo.end(); it++)
		chroms.push_back(it->first);
	if(naturalSort)
		std::sort(chroms.begin(), chroms.end(), naturalChrLess);
	else
		std::sort(chroms.begin(), chroms.end());

	if(verbose)
	{
		time(&tend);
		fprintf(stderr, "Reading has taken %d seconds\n",
		        (int)(difftime(tend, tstart)+0.5));
		cerr << "We have " << nodeCount - 1 << " regions, "
		     << chroms.size() << " chromosomes. Sorting..." << endl;
		time(&tstart);
	}

	// Pass 2 accesses lines in chromosome-sorted order (not file order).
	madvise(mmapBase, mmapSize, MADV_RANDOM);

	// Pass 2: for each chromosome, sort by configured key and print. The
	// default --sort s and --collapse modes don't re-parse — only beg lives
	// in nodes[] and that's all they need. --sort b re-parses end on demand.
	// --sort 5 (needs strand) and --collapse (needs weight) read one extra
	// field per line, and --sort 5 reuses that parse to pull end out for free.
	//
	// At -t > 1 the per-chrom work runs in parallel via std::for_each(par, ...);
	// each chromosome writes to an open_memstream buffer and waits at a
	// producer-consumer barrier for its alphabetical turn before flushing to
	// stdout. --max-mem caps the total bytes a worker can hold mid-flight.
	const bool needStrand = (sortMode == '5') && !fCollapse;

	// Per-chrom worker. Returns 0 on success, 1 on parse error.
	// SortRecSB carries the bucket-sort key (beg, end) plus the line offset so
	// pass 2 can emit without re-reading nodes[]. size_t off bumps the rec to
	// 24 B (16 B + 8 padding to size_t alignment) but is required for
	// files >= 4 GB.
	struct SortRecSB { int beg; int end; size_t off; };
	auto processChrom = [&](int ci, ChromBuf& cb) -> int {
		string2chrInfoT::iterator cit = chrInfo.find(chroms[ci]);

		// Collect node indices for this chromosome (linked-list walk).
		vector<uint32_t> nodeIdx;
		for(uint32_t cur = cit->second.lastRead; cur; cur = nodes[cur].next)
			nodeIdx.push_back(cur);

		if(fCollapse)
		{
			vector<pair<int, float>> wp;
			wp.reserve(nodeIdx.size());
			for(size_t i = 0; i < nodeIdx.size(); i++)
			{
				uint32_t idx = nodeIdx[i];
				const char* linePtr = mmapBase + nodes[idx].off;
				char weight[kWeightBufSize];
				weight[0] = '0'; weight[1] = '\0';
				int dummy_beg = 0, dummy_end = 0;
				char dummy_str = '+';
				const char* chrPtr; int chrLen;
				const char* tailPtr = "";
				int n = parseBedLineFull(linePtr, &chrPtr, &chrLen, &dummy_beg, &dummy_end,
				                         weight, kWeightBufSize, &dummy_str, &tailPtr);
				if(n < 3)
				{
					cerr << "Error in parsing line: " << linePtr << endl;
					return 1;
				}
				char* endptr;
				float w = strtof(weight, &endptr);
				if(endptr == weight)
				{
					cerr << "Malformed weight " << weight << endl;
					return 1;
				}
				wp.emplace_back(nodes[idx].beg, w);
			}
			std::sort(wp.begin(), wp.end(),
			          [](const pair<int, float>& a, const pair<int, float>& b) {
			              return a.first < b.first;
			          });
			size_t i = 0;
			while(i < wp.size())
			{
				int pos = wp[i].first;
				float sum = 0.0f;
				while(i < wp.size() && wp[i].first == pos)
				{
					sum += wp[i].second;
					i++;
				}
				cb.writeFmt("%s\t%d\t%d\t.\t%g\t+\n",
				            chroms[ci].c_str(), pos, pos+1, sum);
			}
			return 0;
		}

		if(needStrand)
		{
			vector<pair<int, uint32_t>> kp;  // (chosenPos, nodeIdx)
			kp.reserve(nodeIdx.size());
			for(size_t i = 0; i < nodeIdx.size(); i++)
			{
				uint32_t idx = nodeIdx[i];
				const char* linePtr = mmapBase + nodes[idx].off;
				char strandChar = '+';
				int dummy_beg = 0, dummy_end = 0;
				char weight[kWeightBufSize]; weight[0] = '0'; weight[1] = '\0';
				const char* chrPtr; int chrLen;
				const char* tailPtr = "";
				int n = parseBedLineFull(linePtr, &chrPtr, &chrLen, &dummy_beg, &dummy_end,
				                         weight, kWeightBufSize, &strandChar, &tailPtr);
				if(n < 5)
				{
					cerr << "Error in parsing line: " << linePtr << endl;
					return 1;
				}
				int pos = (strandChar == '-') ? dummy_end : nodes[idx].beg;
				kp.emplace_back(pos, idx);
			}
			std::sort(kp.begin(), kp.end());
			for(auto& p : kp)
				cb.writeLineFromMmap(mmapBase + nodes[p.second].off);
			return 0;
		}

		// --sort s and --sort b: copy into a contiguous sort-rec for
		// cache-friendly sort, then emit by .off. --sort s leaves .end
		// uninitialised (cmpSB doesn't read it). --sort b re-parses end
		// from the line for the secondary sort key.
		vector<SortRecSB> sr;
		sr.reserve(nodeIdx.size());
		if(sortMode == 'b')
		{
			for(uint32_t idx : nodeIdx)
			{
				const char* linePtr = mmapBase + nodes[idx].off;
				const char* dchrPtr; int dchrLen;
				int dbeg = 0, lineEnd = 0;
				const char* dtailPtr = "";
				if(parseBedLine3(linePtr, &dchrPtr, &dchrLen, &dbeg, &lineEnd, &dtailPtr) < 3)
				{
					cerr << "Error in parsing line: " << linePtr << endl;
					return 1;
				}
				sr.push_back({nodes[idx].beg, lineEnd, nodes[idx].off});
			}
		}
		else
		{
			for(uint32_t idx : nodeIdx)
				sr.push_back({nodes[idx].beg, 0, nodes[idx].off});
		}

		auto cmpSB = [sortMode](const SortRecSB& a, const SortRecSB& b) {
			if(sortMode == 'b')
			{
				if(a.beg != b.beg) return a.beg < b.beg;
				return a.end < b.end;
			}
			return a.beg < b.beg;
		};
		std::sort(sr.begin(), sr.end(), cmpSB);

		for(const auto& r : sr)
			cb.writeLineFromMmap(mmapBase + r.off);
		return 0;
	};

	// Average BED line length (used to pre-size per-chrom output buffers in
	// the parallel path so the buffer doesn't need to double its way up to
	// the chromosome's final size). nodeCount - 1 is the actual line count
	// (index 0 is reserved as the end-of-list sentinel).
	size_t avgLineBytes = (nodeCount > 1)
		? (mmapSize - bodyStart) / (size_t)(nodeCount - 1) + 1
		: 64;

	if(numThreads <= 1)
	{
		// Serial: ChromBuf with flushTo=stdout, 64 KB scratch reused across
		// chromosomes. fwrite_unlocked goes to stdout when the buffer fills,
		// keeping peak RAM bounded.
		ChromBuf cb;
		cb.cap = 64 * 1024;
		cb.buf = (char*) malloc(cb.cap);
		if(!cb.buf) { perror("malloc"); free(nodes); return 1; }
		cb.flushTo = stdout;
		for(size_t ci = 0; ci < chroms.size(); ci++)
		{
			if(processChrom((int)ci, cb) != 0)
			{
				free(cb.buf);
				free(nodes);
				return 1;
			}
			cb.finalize();
		}
		free(cb.buf);
	}
	else
	{
		// Parallel across chromosomes. Each worker writes its chromosome
		// into a pre-sized ChromBuf; producer-consumer barrier flushes
		// in alphabetical order. Optional --max-mem gate caps concurrent
		// per-chromosome buffers (rough estimate: 12 B sort-rec + ~50 B
		// of output per read = ~62 B per read).
		std::vector<int> chromIdxs(chroms.size());
		for(size_t i = 0; i < chroms.size(); i++) chromIdxs[i] = (int)i;

		std::mutex printMtx;
		std::condition_variable printCv;
		int nextChromToPrint = 0;

		std::mutex memMtx;
		std::condition_variable memCv;
		size_t budgetRemaining = maxMemBytes;

		std::for_each(std::execution::par, chromIdxs.begin(), chromIdxs.end(),
			[&](int ci) {
				// Estimate this chrom's peak working set: count its nodes and
				// budget ~62 B per read (12 B sort-rec + ~50 B output buffer).
				string2chrInfoT::iterator cit = chrInfo.find(chroms[ci]);
				size_t chromCount = 0;
				for(uint32_t cur = cit->second.lastRead; cur; cur = nodes[cur].next)
					chromCount++;
				size_t cost = chromCount * 62;
				size_t effCost = (maxMemBytes > 0) ? std::min(cost, maxMemBytes) : 0;

				if(effCost > 0)
				{
					std::unique_lock<std::mutex> lk(memMtx);
					memCv.wait(lk, [&]{ return budgetRemaining >= effCost; });
					budgetRemaining -= effCost;
				}

				ChromBuf cb;
				cb.flushTo = nullptr;
				// Pre-size to ~exactly the expected output size for this
				// chromosome (sort modes preserve line bytes verbatim, so
				// output ≈ chromCount × avgLineBytes). 5% slack absorbs
				// fluctuations; collapse mode may overshoot slightly but
				// the buffer just over-allocates harmlessly.
				cb.cap = std::max((size_t)4096, chromCount * avgLineBytes + chromCount * avgLineBytes / 20);
				cb.buf = (char*) malloc(cb.cap);
				if(!cb.buf)
				{
					cerr << "malloc " << cb.cap << " bytes for chrom "
					     << chroms[ci] << " failed" << endl;
					exit(1);
				}
				int rc = processChrom(ci, cb);

				if(effCost > 0)
				{
					{
						std::lock_guard<std::mutex> lk(memMtx);
						budgetRemaining += effCost;
					}
					memCv.notify_all();
				}

				if(rc != 0)
				{
					free(cb.buf);
					exit(1);
				}

				// Wait for our turn to flush, then write to stdout.
				{
					std::unique_lock<std::mutex> lk(printMtx);
					printCv.wait(lk, [&]{ return nextChromToPrint == ci; });
					fwrite_unlocked(cb.buf, 1, cb.pos, stdout);
					free(cb.buf);
					nextChromToPrint++;
				}
				printCv.notify_all();
			});
	}

	free(nodes);
	if(verbose)
	{
		time(&tend);
		fprintf(stderr, "Sorting has taken %d seconds\n",
		        (int)(difftime(tend, tstart)+0.5));
	}
	return 0;
}

// ============================================================================
// CLASSIC SORT PATH (default)
// ----------------------------------------------------------------------------
// Pipeline:
//   parse: parseLines<UseMmap> (serial / stdin)  OR  parseChunkMmap (parallel
//          via parseMmapDispatch / parseMmapSerial)
//   sort:  sortIndicesDispatch -> sortIndices<SortMode> -> radixSort64 or
//          std::sort with ReadCmp<SortMode>
//   emit:  flat scan over the sorted index, writing each line directly from
//          the mmap buffer (or reconstructing chr\tbeg\tend from the parsed
//          struct + tail string for the stdin/gzip path).
//
// This section also owns the small write helpers (writeUInt, writeBedLine).
// ============================================================================

// Fast integer-to-buffer writer. Returns number of chars written.
static inline int writeUInt(char* buf, int val)
{
	if(val == 0) { buf[0] = '0'; return 1; }
	char tmp[11];
	int n = 0;
	while(val > 0) { tmp[n++] = '0' + (val % 10); val /= 10; }
	for(int i = 0; i < n; i++) buf[i] = tmp[n - 1 - i];
	return n;
}

// Reconstruct and write a BED line: chr\tbeg\tend[tail]\n
// tail includes the leading \t if fields 4+ exist, or is "" for BED3.
static inline void writeBedLine(const char* chr, int beg, int end,
                                const char* tail, FILE* out)
{
	char buf[32];
	fputs_unlocked(chr, out);
	fputc_unlocked('\t', out);
	int n = writeUInt(buf, beg);
	buf[n++] = '\t';
	n += writeUInt(buf + n, end);
	fwrite_unlocked(buf, 1, n, out);
	if(tail[0]) fputs_unlocked(tail, out);
	fputc_unlocked('\n', out);
}


// Comparator that orders read indices by (chrIdx, sort-key). The SortMode
// template parameter selects the secondary key:
//   's' — start coord only
//   'b' — start, then end
//   '5' — 5'-end (end if strand=='-', else beg)
// `if constexpr` collapses each instantiation to a single int compare per arm
// at compile time, giving the same generated code as the old hand-unrolled lambdas.
template<char SortMode>
struct ReadCmp
{
	seqread* reads;
	bool operator()(uint32_t a, uint32_t b) const
	{
		if(reads[a].chrIdx != reads[b].chrIdx) return reads[a].chrIdx < reads[b].chrIdx;
		if constexpr(SortMode == 'b')
		{
			if(reads[a].beg != reads[b].beg) return reads[a].beg < reads[b].beg;
			return reads[a].end < reads[b].end;
		}
		else if constexpr(SortMode == '5')
		{
			int posA = (reads[a].str == '-') ? reads[a].end : reads[a].beg;
			int posB = (reads[b].str == '-') ? reads[b].end : reads[b].beg;
			return posA < posB;
		}
		else  // 's'
		{
			return reads[a].beg < reads[b].beg;
		}
	}
};

// LSD radix sort by 64-bit key, 8-bit digits, 8 passes. Sorts (key, idx)
// pairs so we can reconstruct the sorted index order. Skips passes where
// the byte is constant across all keys (typical for high bytes when
// chrIdx is small and positions fit in 32 bits).
//
// numThreads <= 1: serial path.
// numThreads >  1: thread-local histograms in parallel, serial 2D prefix sum
// over (T threads x 256 buckets), then parallel scatter. Each thread owns a
// disjoint output range so no atomics are needed.
//
// Aux buffers are scratch and freed before return. After the sort, the final
// index order is written back into `order` (the caller's buffer).
static void radixSort64(uint64_t* keys, uint32_t* order, size_t n, int numThreads)
{
	if(n < 2) return;

	uint64_t* keysAux = (uint64_t*) malloc(n * sizeof(uint64_t));
	uint32_t* idxsAux = (uint32_t*) malloc(n * sizeof(uint32_t));
	if(!keysAux || !idxsAux)
	{
		cerr << "Error: out of memory for radix sort scratch (" << n << " entries)" << endl;
		exit(1);
	}

	// Find which bytes vary across all keys; constant bytes can be skipped.
	uint64_t allOr = 0;
	uint64_t allAnd = ~uint64_t(0);
	for(size_t i = 0; i < n; i++) { allOr |= keys[i]; allAnd &= keys[i]; }
	uint64_t varies = allOr ^ allAnd;

	uint64_t* k0 = keys;
	uint64_t* k1 = keysAux;
	uint32_t* i0 = order;
	uint32_t* i1 = idxsAux;
	int passesDone = 0;

	// Decide on chunk count. A few hundred K of work per thread keeps the
	// scheduling overhead amortised; below that the serial path is faster.
	int T = numThreads;
	if(T < 1) T = 1;
	if(n / 200000 < (size_t)T) T = (int)std::max<size_t>(1, n / 200000);

	if(T <= 1)
	{
		// ---- Serial path ----
		for(int pass = 0; pass < 8; pass++)
		{
			if(((varies >> (pass * 8)) & 0xff) == 0) continue;
			const int shift = pass * 8;

			size_t hist[257] = {0};
			for(size_t i = 0; i < n; i++) hist[((k0[i] >> shift) & 0xff) + 1]++;
			for(int b = 1; b < 257; b++) hist[b] += hist[b - 1];

			for(size_t i = 0; i < n; i++)
			{
				size_t pos = hist[(k0[i] >> shift) & 0xff]++;
				k1[pos] = k0[i];
				i1[pos] = i0[i];
			}
			std::swap(k0, k1);
			std::swap(i0, i1);
			passesDone++;
		}
	}
	else
	{
		// ---- Parallel path ----
		std::vector<size_t> chunkStart(T + 1);
		for(int t = 0; t <= T; t++)
			chunkStart[t] = n * (size_t)t / (size_t)T;

		std::vector<int> tIds(T);
		for(int t = 0; t < T; t++) tIds[t] = t;

		// Per-thread histograms: hist[t * 256 + b].
		std::vector<size_t> hist((size_t)T * 256);
		// Per-thread, per-bucket write cursors.
		std::vector<size_t> cursor((size_t)T * 256);

		for(int pass = 0; pass < 8; pass++)
		{
			if(((varies >> (pass * 8)) & 0xff) == 0) continue;
			const int shift = pass * 8;

			// Phase 1: parallel local histograms.
			std::for_each(std::execution::par, tIds.begin(), tIds.end(), [&](int t) {
				size_t* h = hist.data() + (size_t)t * 256;
				for(int b = 0; b < 256; b++) h[b] = 0;
				size_t e = chunkStart[t + 1];
				for(size_t i = chunkStart[t]; i < e; i++)
					h[(k0[i] >> shift) & 0xff]++;
			});

			// Phase 2: serial 2D prefix sum (T*256 cells = ~6 KB at T=22).
			// cursor[t * 256 + b] = starting write position for thread t's
			// bucket b in the output buffer.
			size_t base = 0;
			for(int b = 0; b < 256; b++)
			{
				for(int t = 0; t < T; t++)
				{
					cursor[(size_t)t * 256 + b] = base;
					base += hist[(size_t)t * 256 + b];
				}
			}

			// Phase 3: parallel scatter. Each thread's writes go into a disjoint
			// per-(thread, bucket) range, so no contention.
			std::for_each(std::execution::par, tIds.begin(), tIds.end(), [&](int t) {
				size_t* c = cursor.data() + (size_t)t * 256;
				size_t e = chunkStart[t + 1];
				for(size_t i = chunkStart[t]; i < e; i++)
				{
					size_t pos = c[(k0[i] >> shift) & 0xff]++;
					k1[pos] = k0[i];
					i1[pos] = i0[i];
				}
			});

			std::swap(k0, k1);
			std::swap(i0, i1);
			passesDone++;
		}
	}

	// If we did an odd number of passes, the final data is in the aux buffer;
	// copy it back so the caller's `order` holds the result.
	if(passesDone & 1)
		memcpy(order, i0, n * sizeof(uint32_t));

	free(keysAux);
	free(idxsAux);
}

// Build 64-bit radix keys for SortMode 's' or '5'.
// Layout: high 32 bits = chrIdx, low 32 bits = position (beg or 5'-end).
template<char SortMode>
static void buildRadixKeys(uint64_t* keys, const uint32_t* order, size_t n, const seqread* reads)
{
	static_assert(SortMode == 's' || SortMode == '5',
	              "radix path only supports SortMode 's' or '5'");
	for(size_t i = 0; i < n; i++)
	{
		uint32_t idx = order[i];
		uint32_t pos;
		if constexpr(SortMode == '5')
			pos = (reads[idx].str == '-') ? (uint32_t)reads[idx].end : (uint32_t)reads[idx].beg;
		else
			pos = (uint32_t)reads[idx].beg;
		keys[i] = ((uint64_t)reads[idx].chrIdx << 32) | pos;
	}
}

// Radix path thresholds. Single-thread: ~100k is where radix amortises its
// histogram+double-buffer cost vs std::sort. Multi-thread: parallel radix
// has more setup (T*256 histograms, two parallel_for syncs per pass), so
// stay on the heavily-tuned std::sort(par) below ~3M reads.
static constexpr size_t RADIX_SORT_THRESHOLD       = 100000;
static constexpr size_t RADIX_SORT_THRESHOLD_PAR   = 3000000;

template<char SortMode>
static void sortIndices(uint32_t* order, size_t n, seqread* reads, int numThreads)
{
	// Radix path: large enough, sortMode supports a single 32-bit position key.
	// --sort b's three keys (chrIdx, beg, end) don't fit in 64 bits without
	// lossy squeezing, so it never reaches this template instantiation — the
	// dispatcher in main() routes --sort b to sortIndicesPerChromB instead.
	// Same key-packing works for both single-thread and multi-thread;
	// radixSort64 internally dispatches to a serial or parallel histogram +
	// scatter pipeline.
	if constexpr(SortMode == 's' || SortMode == '5')
	{
		const size_t threshold = (numThreads > 1) ? RADIX_SORT_THRESHOLD_PAR : RADIX_SORT_THRESHOLD;
		if(n >= threshold)
		{
			uint64_t* keys = (uint64_t*) malloc(n * sizeof(uint64_t));
			if(!keys)
			{
				cerr << "Error: out of memory for radix keys (" << n << " entries)" << endl;
				exit(1);
			}
			buildRadixKeys<SortMode>(keys, order, n, reads);
			radixSort64(keys, order, n, numThreads);
			free(keys);
			return;
		}
	}

	// Fallback: comparator-based std::sort. Used for --sort b at any thread
	// count, and for small n where radix's scheduling overhead would dominate.
	ReadCmp<SortMode> cmp{reads};
	if(numThreads == 1)
		std::sort(order, order + n, cmp);
	else
		std::sort(std::execution::par, order, order + n, cmp);
}

// Dispatch the templated sort by the runtime sortMode char.
static void sortIndicesDispatch(uint32_t* order, size_t n, seqread* reads, char sortMode, int numThreads)
{
	switch(sortMode)
	{
		case 'b': sortIndices<'b'>(order, n, reads, numThreads); break;
		case '5': sortIndices<'5'>(order, n, reads, numThreads); break;
		default:  sortIndices<'s'>(order, n, reads, numThreads); break;
	}
}

// Per-chromosome radix sort for --sort b. The classic-path order[] array is
// already laid out in chromosome order (caller populates one chrom at a
// time), so we can sort each chrom's slice independently. Within a single
// chromosome the key is (beg << 32) | end — fits in 64 bits cleanly, no
// chrIdx packing needed. Below RADIX_SORT_THRESHOLD per chrom we fall
// through to std::sort with the comparator.
//
// Chromosomes are processed in parallel via std::for_each(par) at -t > 1;
// each worker runs serial radix on its assigned chromosome. This avoids
// the per-pass synchronisation cost of parallel-radix-within-one-chrom and
// scales to all available worker threads when there are many chromosomes.
static void sortIndicesPerChromB(uint32_t* order,
                                 const std::vector<size_t>& chromStart,
                                 seqread* reads, int numThreads)
{
	size_t numChroms = chromStart.size() - 1;

	auto sortChrom = [&](size_t ci) {
		size_t start = chromStart[ci];
		size_t n = chromStart[ci + 1] - start;
		if(n < 2) return;

		if(n < RADIX_SORT_THRESHOLD)
		{
			ReadCmp<'b'> cmp{reads};
			std::sort(order + start, order + start + n, cmp);
			return;
		}

		// Build 64-bit keys (beg, end) for this chrom's slice.
		uint64_t* keys = (uint64_t*) malloc(n * sizeof(uint64_t));
		if(!keys)
		{
			cerr << "Error: out of memory for radix keys (" << n << " entries)" << endl;
			exit(1);
		}
		for(size_t i = 0; i < n; i++)
		{
			uint32_t idx = order[start + i];
			keys[i] = ((uint64_t)(uint32_t)reads[idx].beg << 32) |
			          (uint64_t)(uint32_t)reads[idx].end;
		}
		// Serial radix per chrom (parallelism is across chroms, not within).
		radixSort64(keys, order + start, n, 1);
		free(keys);
	};

	if(numThreads <= 1 || numChroms < 2)
	{
		for(size_t ci = 0; ci < numChroms; ci++) sortChrom(ci);
	}
	else
	{
		std::vector<size_t> chromIdxs(numChroms);
		for(size_t i = 0; i < numChroms; i++) chromIdxs[i] = i;
		std::for_each(std::execution::par, chromIdxs.begin(), chromIdxs.end(), sortChrom);
	}
}

// Parsing loop, templated on UseMmap to eliminate per-line branch.
// Stores only the "tail" (fields 4+) in reads[].line (mmap: full line pointer;
// stdin/gzip: only tail copied to arena, saving ~50% arena memory).
template<bool UseMmap>
static void parseLines(
    char* mmapBase, size_t mmapSize,
    FILE* fh,
    seqread*& reads, size_t& readCount, size_t& currMaxReads,
    string2chrInfoT& chrInfo, string2chrInfoT::iterator& thisChrIt,
    Arena* arena,
    int fCollapse, char sortMode)
{
	char* mmapCur = mmapBase;
	char* mmapLim = mmapBase + mmapSize;
	// getline() reuses and grows this buffer across calls, so there's no
	// fixed line-length limit on the stdin/gzip path. Only allocated when
	// UseMmap is false; freed at end of function.
	char* lineBuf = NULL;
	size_t lineBufCap = 0;
	const bool needExtra = fCollapse || (sortMode == '5');

	while(1)
	{
		char* linePtr;

		if(UseMmap)
		{
			if(mmapCur >= mmapLim) break;
			linePtr = mmapCur;
			char* nl = (char*) memchr(mmapCur, '\n', (size_t)(mmapLim - mmapCur));
			if(nl)
			{
				if(nl > mmapCur && *(nl - 1) == '\r') *(nl - 1) = '\0';
				*nl = '\0';
				mmapCur = nl + 1;
			}
			else
			{
				mmapCur = mmapLim;
			}
		}
		else
		{
			ssize_t glen = getline(&lineBuf, &lineBufCap, fh);
			if(glen == -1) break;  // EOF or error
			linePtr = lineBuf;
			// Strip trailing \r\n or \n (getline keeps the newline if present).
			if(glen >= 2 && lineBuf[glen-2] == '\r' && lineBuf[glen-1] == '\n')
				lineBuf[glen-2] = '\0';
			else if(glen >= 1 && lineBuf[glen-1] == '\n')
				lineBuf[glen-1] = '\0';
		}

		// Pass through BED header lines (track/browser/# comments) directly.
		if(linePtr[0] == '#' ||
		   strncmp(linePtr, "track ", 6) == 0 ||
		   strncmp(linePtr, "browser ", 8) == 0)
		{
			fputs_unlocked(linePtr, stdout);
			fputc_unlocked('\n', stdout);
			continue;
		}

		if(readCount == currMaxReads)
		{
			currMaxReads *= 2;
			seqread* tmp = (seqread*) realloc(reads, currMaxReads * sizeof(seqread));
			if(!tmp) { fprintf(stderr, "Error: out of memory\n"); exit(1); }
			reads = tmp;
		}
		// reads[].next, the order[] array, and radix-sort internals are all
		// uint32_t. Index 0 is reserved as the end-of-list sentinel, so
		// usable range is 1..UINT32_MAX-1.
		if(readCount >= (size_t)UINT32_MAX)
		{
			fprintf(stderr, "Error: %zu reads exceeds the 4.29 B limit "
			        "imposed by the 32-bit per-read index. Split your input.\n",
			        readCount);
			exit(1);
		}

		int beg = 0;
		int end = 0;
		char strandChar = '+';
		char weight[kWeightBufSize];
		weight[0] = '0'; weight[1] = '\0';
		const char* chrPtr;
		int chrLen;
		const char* tailPtr = "";
		int numArgsRead;

		if(needExtra)
		{
			numArgsRead = parseBedLineFull(linePtr, &chrPtr, &chrLen,
			                              &beg, &end, weight, kWeightBufSize,
			                              &strandChar, &tailPtr);
		}
		else
		{
			numArgsRead = parseBedLine3(linePtr, &chrPtr, &chrLen, &beg, &end, &tailPtr);
		}

		if(numArgsRead >= 3)
		{
			if(beg < 0)
			{
				cerr << "Error: negative coordinates in the bed file\n" << linePtr << endl;
				exit(1);
			}
			if(!fCollapse)
			{
				if(UseMmap)
				{
					reads[readCount].line = linePtr;
				}
				else
				{
					// stdin/gzip: store only tail to save arena memory
					size_t tlen = strlen(tailPtr);
					reads[readCount].line = arena->alloc(tailPtr, tlen + 1);
				}
			}
			else
			{
				size_t wlen = strlen(weight);
				reads[readCount].line = arena->alloc(weight, wlen + 1);
			}
			reads[readCount].beg = beg;
			reads[readCount].end = end;
			reads[readCount].str = strandChar;

			// Fast same-chr check: compare length first, then bytes.
			if((int)thisChrIt->first.size() != chrLen ||
			   memcmp(thisChrIt->first.data(), chrPtr, chrLen) != 0)
			{
				// tryEmplaceByPtr avoids std::string construction on the common
				// "chrom already exists" path; only allocates on a true insert.
				auto ins = chrInfo.tryEmplaceByPtr(chrPtr, chrLen);
				thisChrIt = ins.first;
				if(ins.second)
					reads[readCount].next = 0;
				else
					reads[readCount].next = thisChrIt->second.lastRead;
			}
			else
			{
				reads[readCount].next = thisChrIt->second.lastRead;
			}
			thisChrIt->second.lastRead = (uint32_t)readCount;
		}
		else
		{
			cerr << "Error in parsing line: " << linePtr << endl
				<< "Perhaps this line is malformed?" << endl;
			exit(1);
		}
		readCount ++;
	}
	free(lineBuf);
}

// ============================================================================
// Parallel mmap parsing (used when --threads > 1 on a non-gzip mmap input)
// ----------------------------------------------------------------------------
// Mmap is split into N newline-aligned chunks. A first parallel pass counts
// newlines per chunk so we can size the global `reads` array exactly;
// prefix-sums then give each chunk a base index. A second parallel pass
// parses each chunk into its slot, building per-chunk per-chrom partial
// linked lists. A serial merge step concatenates them.
//
// Header lines (track/browser/#) are emitted to stdout *before* chunking,
// in a leading pre-pass over the mmap (only header lines at the start of
// the file are recognized; per BED convention headers are leading).
// ============================================================================

struct ChunkChrPartial
{
	uint32_t head;   // global idx of head (most recently inserted read)
	uint32_t tail;   // global idx of tail (read whose .next is 0)
};

typedef ChrNameMap<ChunkChrPartial> ChunkChrMap;

struct ChunkResult
{
	uint32_t firstIdx;
	uint32_t count;
	ChunkChrMap chr;
};

// Parse [start, end) of the mmap into reads[base+1 .. base+count].
// All indices stored (head/tail and reads[i].next) are global, so no
// rebase pass is needed. Mutates the chunk in place: replaces '\n' (and
// preceding '\r') with '\0' to NUL-terminate each line.
//
// If fCollapse is set, the per-chunk `arena` is used to copy weight strings
// out of the (possibly clobbered) line buffer; reads[i].line then points at
// the arena copy. arena must be non-null whenever fCollapse is set.
static void parseChunkMmap(char* start, char* end,
                           seqread* reads, uint32_t base,
                           ChunkResult& result,
                           int fCollapse, char sortMode,
                           Arena* arena)
{
	const bool needExtra = fCollapse || (sortMode == '5');
	char* mmapCur = start;
	char* mmapLim = end;
	uint32_t idx = base + 1;
	auto thisChrIt = result.chr.end();

	while(mmapCur < mmapLim)
	{
		char* linePtr = mmapCur;
		char* nl = (char*) memchr(mmapCur, '\n', (size_t)(mmapLim - mmapCur));
		if(nl)
		{
			if(nl > mmapCur && *(nl - 1) == '\r') *(nl - 1) = '\0';
			*nl = '\0';
			mmapCur = nl + 1;
		}
		else
		{
			mmapCur = mmapLim;
		}

		int beg = 0, lineEnd = 0;
		char strandChar = '+';
		char weight[kWeightBufSize];
		weight[0] = '0'; weight[1] = '\0';
		const char* chrPtr;
		int chrLen;
		const char* tailPtr = "";
		int numArgsRead;

		if(needExtra)
		{
			numArgsRead = parseBedLineFull(linePtr, &chrPtr, &chrLen,
			                               &beg, &lineEnd, weight, kWeightBufSize,
			                               &strandChar, &tailPtr);
		}
		else
		{
			numArgsRead = parseBedLine3(linePtr, &chrPtr, &chrLen, &beg, &lineEnd, &tailPtr);
		}

		if(numArgsRead < 3)
		{
			cerr << "Error in parsing line: " << linePtr << endl
				 << "Perhaps this line is malformed?" << endl;
			exit(1);
		}
		if(beg < 0)
		{
			cerr << "Error: negative coordinates in the bed file\n" << linePtr << endl;
			exit(1);
		}

		// Collapse mode stores only the weight (used in the output sum); the line
		// pointer goes into a per-chunk arena so we don't depend on the mmap
		// buffer beyond the parser. Non-collapse: zero-copy mmap pointer.
		if(fCollapse)
		{
			size_t wlen = strlen(weight);
			reads[idx].line = arena->alloc(weight, wlen + 1);
		}
		else
		{
			reads[idx].line = linePtr;
		}
		reads[idx].beg = beg;
		reads[idx].end = lineEnd;
		reads[idx].str = strandChar;

		// Same-chr fast path within this chunk.
		if(thisChrIt != result.chr.end() &&
		   (int)thisChrIt->first.size() == chrLen &&
		   memcmp(thisChrIt->first.data(), chrPtr, chrLen) == 0)
		{
			reads[idx].next = thisChrIt->second.head;
			thisChrIt->second.head = idx;
		}
		else
		{
			auto ins = result.chr.tryEmplaceByPtr(chrPtr, chrLen);
			thisChrIt = ins.first;
			if(ins.second)
			{
				reads[idx].next = 0;
				thisChrIt->second.head = idx;
				thisChrIt->second.tail = idx;
			}
			else
			{
				reads[idx].next = thisChrIt->second.head;
				thisChrIt->second.head = idx;
			}
		}
		idx++;
	}

	result.firstIdx = base + 1;
	result.count = idx - (base + 1);
}

// Single-thread mmap parser. Wraps parseLines<true> with the same return shape
// as parseMmapParallel so the dispatch site can ignore which one ran. Skips
// the pre-count pass entirely (the legacy parser estimates capacity from the
// first 64 KB and grows via realloc — one less full pass over the mmap).
//
// Body parsing starts at `bodyStart` (after any leading header lines that
// the caller has already emitted to stdout).
static seqread* parseMmapSerial(char* mmapBase, size_t bodyStart, size_t mmapSize,
                                int fCollapse, char sortMode,
                                Arena* arena,
                                string2chrInfoT& chrInfo,
                                size_t& outReadCount)
{
	size_t currMaxReads = 1024;
	size_t bodySize = mmapSize - bodyStart;
	size_t scanSize = bodySize < 65536 ? bodySize : 65536;
	size_t nLines = 0;
	for(size_t i = 0; i < scanSize; i++)
		if(mmapBase[bodyStart + i] == '\n') nLines++;
	if(nLines > 0)
		currMaxReads = (size_t)((double)bodySize / (double)scanSize * (double)nLines * 1.1);
	if(currMaxReads < 1024) currMaxReads = 1024;

	seqread* reads = (seqread*) malloc(currMaxReads * sizeof(seqread));
	if(!reads)
	{
		cerr << "Error: out of memory allocating reads array" << endl;
		exit(1);
	}

	size_t readCount = 1;  // slot 0 is the end-of-list sentinel
	auto ins = chrInfo.insert(std::make_pair(std::string(""), chrInfoT()));
	auto thisChrIt = ins.first;
	parseLines<true>(mmapBase + bodyStart, bodySize, NULL,
	                 reads, readCount, currMaxReads,
	                 chrInfo, thisChrIt, arena,
	                 fCollapse, sortMode);
	chrInfo.erase("");
	outReadCount = readCount - 1;
	return reads;
}

// Top-level mmap parser. Emits any leading header lines to stdout, then
// either runs parseMmapSerial (N==1) or chunks the body and parses in
// parallel. Returns the malloc'd reads array; caller frees.
static seqread* parseMmapDispatch(char* mmapBase, size_t mmapSize,
                                  int fCollapse, char sortMode,
                                  int numThreads,
                                  Arena* arena,
                                  std::vector<Arena*>* chunkArenas,
                                  string2chrInfoT& chrInfo,
                                  size_t& outReadCount)
{
	// Pre-pass: emit leading header lines (#, track, browser) to stdout and
	// advance bodyStart past them. BED convention puts headers at the top.
	size_t bodyStart = 0;
	while(bodyStart < mmapSize)
	{
		char* lineStart = mmapBase + bodyStart;
		char* nl = (char*) memchr(lineStart, '\n', mmapSize - bodyStart);
		size_t lineLen = (size_t)((nl ? nl : mmapBase + mmapSize) - lineStart);
		bool isHeader = (lineLen > 0 && lineStart[0] == '#') ||
		                (lineLen >= 6 && memcmp(lineStart, "track ",   6) == 0) ||
		                (lineLen >= 8 && memcmp(lineStart, "browser ", 8) == 0);
		if(!isHeader) break;
		fwrite_unlocked(lineStart, 1, lineLen, stdout);
		fputc_unlocked('\n', stdout);
		bodyStart = (nl ? (size_t)(nl - mmapBase) + 1 : mmapSize);
	}

	int N = numThreads;
	if(N < 1) N = 1;
	{
		size_t bytesPerChunk = 256 * 1024;
		size_t bodySize = mmapSize - bodyStart;
		int byBytes = (int)((bodySize + bytesPerChunk - 1) / bytesPerChunk);
		if(byBytes < 1) byBytes = 1;
		if(N > byBytes) N = byBytes;
	}

	// Single-thread fast path uses the legacy realloc-grow parser to avoid
	// a separate pre-count pass over the mmap.
	if(N == 1)
		return parseMmapSerial(mmapBase, bodyStart, mmapSize,
		                       fCollapse, sortMode, arena,
		                       chrInfo, outReadCount);

	// In --collapse mode each chunk needs its own Arena to hold weight strings
	// (mutex-sharing one would serialise the parse). The arenas are kept alive
	// by the caller's `chunkArenas` vector for the duration of the program.
	if(fCollapse)
	{
		chunkArenas->resize(N);
		for(int i = 0; i < N; i++)
			(*chunkArenas)[i] = new Arena(4UL * 1024 * 1024, 64UL * 1024 * 1024);
	}

	// Newline-aligned chunk boundaries within [bodyStart, mmapSize).
	std::vector<size_t> chunkStart(N + 1);
	chunkStart[0] = bodyStart;
	chunkStart[N] = mmapSize;
	for(int i = 1; i < N; i++)
	{
		size_t guess = bodyStart + ((mmapSize - bodyStart) / (size_t)N) * (size_t)i;
		while(guess < mmapSize && mmapBase[guess] != '\n') guess++;
		if(guess < mmapSize) guess++;
		chunkStart[i] = guess;
	}

	// Pass 1: count newlines per chunk in parallel.
	std::vector<size_t> chunkCount(N, 0);
	std::vector<int> indices(N);
	for(int i = 0; i < N; i++) indices[i] = i;
	std::for_each(std::execution::par, indices.begin(), indices.end(),
		[&](int i) {
			size_t a = chunkStart[i], b = chunkStart[i + 1];
			size_t c = 0;
			const char* p = mmapBase + a;
			const char* pe = mmapBase + b;
			while(p < pe)
			{
				const char* q = (const char*) memchr(p, '\n', (size_t)(pe - p));
				if(!q) break;
				c++;
				p = q + 1;
			}
			chunkCount[i] = c;
		});
	// Trailing-no-newline: the last chunk has one extra line.
	if(mmapSize > bodyStart && mmapBase[mmapSize - 1] != '\n')
		chunkCount[N - 1]++;

	std::vector<size_t> chunkBase(N + 1);
	chunkBase[0] = 0;
	for(int i = 0; i < N; i++) chunkBase[i + 1] = chunkBase[i] + chunkCount[i];
	size_t totalReads = chunkBase[N];

	// reads[].next, the order[] array, and radix-sort internals are all
	// uint32_t. Index 0 is reserved as the end-of-list sentinel, so usable
	// range is 1..UINT32_MAX-1.
	if(totalReads >= (size_t)UINT32_MAX)
	{
		cerr << "Error: " << totalReads << " reads exceeds the 4.29 B limit "
		     << "imposed by the 32-bit per-read index. Split your input." << endl;
		exit(1);
	}

	seqread* reads = (seqread*) malloc((totalReads + 1) * sizeof(seqread));
	if(!reads)
	{
		cerr << "Error: out of memory allocating reads array (" << totalReads << " entries)" << endl;
		exit(1);
	}

	std::vector<ChunkResult> chunkRes(N);
	std::for_each(std::execution::par, indices.begin(), indices.end(),
		[&](int i) {
			Arena* a = fCollapse ? (*chunkArenas)[i] : NULL;
			parseChunkMmap(mmapBase + chunkStart[i], mmapBase + chunkStart[i + 1],
			               reads, chunkBase[i], chunkRes[i],
			               fCollapse, sortMode, a);
		});

	for(int i = 0; i < N; i++)
	{
		if(chunkRes[i].count != chunkCount[i])
		{
			cerr << "Internal error: chunk " << i << " parsed " << chunkRes[i].count
			     << " reads, expected " << chunkCount[i] << endl;
			exit(1);
		}
	}

	// Merge: concatenate per-chunk per-chrom lists into the global chrInfo.
	for(int t = 0; t < N; t++)
	{
		for(auto& kv : chunkRes[t].chr)
		{
			const std::string& name = kv.first;
			const ChunkChrPartial& partial = kv.second;
			auto ins = chrInfo.try_emplace(name);
			chrInfoT& info = ins.first->second;
			if(ins.second)
			{
				info.lastRead = partial.head;
			}
			else
			{
				reads[partial.tail].next = info.lastRead;
				info.lastRead = partial.head;
			}
		}
	}

	outReadCount = totalReads;
	return reads;
}

#ifdef WITH_BAM
// ============================================================================
// BAM PATH (--bam / *.bam input, optional, opt-in at build time via WITH_BAM)
// ----------------------------------------------------------------------------
// In-RAM coordinate sort using htslib. Reads every record into a vector,
// builds a (tid, pos) packed-uint64 key per record, runs the existing
// radixSort64 for the index permutation, then writes the records out as BAM
// in sorted order via sam_write1. hts_set_threads enables htslib's worker
// pool for parallel BGZF compression/decompression.
//
// This is the "Approach A" MVP: simple, beats samtools on RAM-resident BAMs
// (parallel radix vs. ksort comparator), but holds all records in memory so
// it tops out around half-RAM input size. Larger BAMs need spill+merge
// (samtools-clone) which is intentionally out of scope here.
// ============================================================================
static int bamSortAndEmit(const std::string& inputFile,
                          const std::string& outputFile,
                          int numThreads,
                          bool verbose)
{
	samFile* in = sam_open(inputFile.c_str(), "r");
	if(!in)
	{
		std::cerr << "Error: cannot open BAM " << inputFile << std::endl;
		return 1;
	}
	if(numThreads > 1) hts_set_threads(in, numThreads);

	sam_hdr_t* hdr = sam_hdr_read(in);
	if(!hdr)
	{
		std::cerr << "Error: cannot read BAM header from " << inputFile << std::endl;
		sam_close(in);
		return 1;
	}

	time_t tstart, tend;
	if(verbose) time(&tstart);

	std::vector<bam1_t*> records;
	while(true)
	{
		bam1_t* b = bam_init1();
		int r = sam_read1(in, hdr, b);
		if(r < -1)
		{
			std::cerr << "Error: truncated/corrupt BAM at record "
			          << records.size() << std::endl;
			bam_destroy1(b);
			sam_hdr_destroy(hdr);
			sam_close(in);
			return 1;
		}
		if(r == -1) { bam_destroy1(b); break; }
		records.push_back(b);
	}
	sam_close(in);

	if(verbose)
	{
		time(&tend);
		std::cerr << "BAM read: " << records.size() << " records in "
		          << (long)(tend - tstart) << " s" << std::endl;
		time(&tstart);
	}

	// Update header SO before write; works for empty input too.
	sam_hdr_update_hd(hdr, "SO", "coordinate");

	const char* outPath = outputFile.empty() ? "-" : outputFile.c_str();
	samFile* out = sam_open(outPath, "wb");
	if(!out)
	{
		std::cerr << "Error: cannot open output BAM " << outPath << std::endl;
		sam_hdr_destroy(hdr);
		for(auto* b : records) bam_destroy1(b);
		return 1;
	}
	if(numThreads > 1) hts_set_threads(out, numThreads);

	size_t n = records.size();
	uint32_t* order = NULL;
	if(n > 0)
	{
		uint64_t* keys = (uint64_t*) malloc(n * sizeof(uint64_t));
		order = (uint32_t*) malloc(n * sizeof(uint32_t));
		if(!keys || !order)
		{
			std::cerr << "Error: out of memory for BAM sort index ("
			          << n << " entries)" << std::endl;
			sam_close(out);
			sam_hdr_destroy(hdr);
			for(auto* b : records) bam_destroy1(b);
			return 1;
		}
		for(size_t i = 0; i < n; i++)
		{
			// tid == -1 (unmapped) wraps to UINT32_MAX, which sorts after every
			// real reference id. pos is 0-based int64_t in modern htslib; cast
			// to uint32_t — fine for any standard reference (max chrom < 2^31).
			uint32_t tid = (uint32_t) records[i]->core.tid;
			uint32_t pos = (uint32_t) records[i]->core.pos;
			keys[i] = ((uint64_t)tid << 32) | pos;
			order[i] = (uint32_t) i;
		}
		radixSort64(keys, order, n, numThreads);
		free(keys);

		if(verbose)
		{
			time(&tend);
			std::cerr << "BAM sort: " << (long)(tend - tstart) << " s" << std::endl;
			time(&tstart);
		}
	}

	if(sam_hdr_write(out, hdr) < 0)
	{
		std::cerr << "Error: cannot write BAM header" << std::endl;
		sam_close(out);
		sam_hdr_destroy(hdr);
		free(order);
		for(auto* b : records) bam_destroy1(b);
		return 1;
	}
	for(size_t i = 0; i < n; i++)
	{
		if(sam_write1(out, hdr, records[order[i]]) < 0)
		{
			std::cerr << "Error: cannot write BAM record " << i << std::endl;
			sam_close(out);
			sam_hdr_destroy(hdr);
			free(order);
			for(auto* b : records) bam_destroy1(b);
			return 1;
		}
	}
	sam_close(out);

	if(verbose && n > 0)
	{
		time(&tend);
		std::cerr << "BAM write: " << (long)(tend - tstart) << " s" << std::endl;
	}

	for(auto* b : records) bam_destroy1(b);
	free(order);
	sam_hdr_destroy(hdr);
	return 0;
}
#endif // WITH_BAM

// ============================================================================
// EXTERNAL MERGE SORT PATH (--external-merge)
// ----------------------------------------------------------------------------
// Streaming sort for inputs > RAM. Pass 1 reads the file in chunks, fills an
// in-RAM run buffer up to a memory budget, sorts with radixSort64, and writes
// a binary run file to a temp directory. Pass 2 (merge) opens all run files,
// maintains a min-heap of (key, runIdx), pops + emits + refills.
//
// Run-file format (codec=0=raw — compression codecs added in later commits):
//   Header:
//     magic "PSBR" (4 B), version u8 (=1), codec u8, flags u16,
//     chrCount u32, chrs: [u16 nameLen | nameLen bytes] × chrCount
//   Blocks (sequential, one or more):
//     u32 uncompLen, u32 compLen, u32 numRecords, payload (compLen bytes)
//   Footer:
//     magic "PSBE" (4 B), u32 blockCount, u64 totalRecords, u32 padding
// Each block payload is a sequence of records:
//     u16 chrIdx, i32 beg, i32 end, u8 strand, u16 tailLen, tailLen bytes tail
// Tail = BED fields 4+ (no chr name; chrIdx indexes into the run's chr dict).
// ============================================================================

enum ExtCodec : uint8_t {
	EXT_RAW  = 0,
	EXT_LZ4  = 1,
	EXT_ZSTD = 2,
	EXT_RANS0 = 3,
	EXT_RANS1 = 4,
};

static const char* extCodecName(ExtCodec c) {
	switch(c) {
		case EXT_RAW:   return "raw";
		case EXT_LZ4:   return "lz4";
		case EXT_ZSTD:  return "zstd";
		case EXT_RANS0: return "rans0";
		case EXT_RANS1: return "rans1";
	}
	return "?";
}

static int parseExtCodec(const std::string& s, ExtCodec* out) {
	if(s == "raw")   { *out = EXT_RAW;   return 0; }
	if(s == "lz4")   { *out = EXT_LZ4;   return 0; }
	if(s == "zstd")  { *out = EXT_ZSTD;  return 0; }
#ifdef WITH_RANS
	if(s == "rans0") { *out = EXT_RANS0; return 0; }
	if(s == "rans1") { *out = EXT_RANS1; return 0; }
#else
	if(s == "rans0" || s == "rans1") {
		std::cerr << "Error: codec '" << s << "' requires a WITH_BAM=1 build "
		          << "(htscodecs rANS is bundled into libhts.a)" << std::endl;
		return -1;
	}
#endif
	return -1;
}

// 1 MiB uncompressed block target — chosen so parallel zstd
// (ZSTD_c_nbWorkers) has enough payload per block to engage workers.
// Smaller blocks (e.g., 64 KB) silently degrade to single-stream zstd.
// Bigger blocks marginally improve ratio but cost RAM in the writer.
constexpr size_t kExtBlockTargetBytes = 1024 * 1024;

// One parsed record in the in-RAM run buffer. Tail bytes live in a separate
// Arena (allocated via Arena::alloc which returns a pointer); the entry just
// holds the pointer + length so the array stays cheap to sort.
struct ExtEntry {
	uint16_t    chrIdx;
	int32_t     beg;
	int32_t     end;
	uint8_t     strand;
	uint16_t    tailLen;
	const char* tailPtr;
};

// Compress one block. Returns compressed size; writes to outBuf which must
// have at least extCompressBound(codec, inLen) bytes available.
//
// numThreads > 1: enables parallel zstd via ZSTD_c_nbWorkers (other codecs
// have no built-in MT in our setup and ignore it). Block size needs to be
// >= ~nbWorkers x 128 KB for zstd MT to actually engage; below that zstd
// silently falls back to single-stream. The writer's block target
// kExtBlockTargetBytes = 1 MiB makes zstd MT useful for nbWorkers <= 8.
static size_t extCompressBlock(ExtCodec codec, const uint8_t* in, size_t inLen,
                               uint8_t* outBuf, size_t maxOut, int numThreads) {
	switch(codec) {
		case EXT_RAW:
			memcpy(outBuf, in, inLen);
			return inLen;
		case EXT_LZ4: {
			int n = LZ4_compress_default((const char*)in, (char*)outBuf,
			                             (int)inLen, (int)maxOut);
			if(n <= 0) {
				std::cerr << "Error: LZ4_compress_default failed (in="
				          << inLen << ", maxOut=" << maxOut << ")" << std::endl;
				exit(1);
			}
			return (size_t)n;
		}
		case EXT_ZSTD: {
			// Empirically, ZSTD_c_nbWorkers > 0 *regresses* on our 1 MiB
			// blocks (200 M / 4 G bench: extmerge -t 8 = 135 s vs -t 1
			// = 117 s with MT enabled; turning it off lets the radix-sort
			// MT win cleanly). zstd's MT mode chunks the input into
			// ~128 KB pieces per worker, but the dispatch+sync cost at
			// our block size exceeds the parallelism benefit. Bigger
			// blocks (>= 16 MiB) would help, at the cost of writer RAM
			// and emit granularity. Keep zstd single-stream for now.
			(void)numThreads;
			size_t n = ZSTD_compress(outBuf, maxOut, in, inLen, 1);
			if(ZSTD_isError(n)) {
				std::cerr << "Error: ZSTD_compress failed: "
				          << ZSTD_getErrorName(n) << std::endl;
				exit(1);
			}
			return n;
		}
#ifdef WITH_RANS
		case EXT_RANS0:
		case EXT_RANS1: {
			int order = (codec == EXT_RANS0) ? 0 : 1;
			unsigned int outSize = (unsigned int)maxOut;
			unsigned char* r = rans_compress_to_4x16(
				(unsigned char*)in, (unsigned int)inLen,
				(unsigned char*)outBuf, &outSize, order);
			if(!r) {
				std::cerr << "Error: rans_compress_to_4x16 (order " << order
				          << ") failed" << std::endl;
				exit(1);
			}
			return (size_t)outSize;
		}
#endif
		default:
			std::cerr << "Error: codec " << extCodecName(codec)
			          << " not available in this build" << std::endl;
			exit(1);
	}
}

static void extDecompressBlock(ExtCodec codec, const uint8_t* in, size_t compLen,
                               uint8_t* out, size_t uncompLen) {
	switch(codec) {
		case EXT_RAW:
			if(compLen != uncompLen) {
				std::cerr << "Error: raw block size mismatch" << std::endl;
				exit(1);
			}
			memcpy(out, in, uncompLen);
			return;
		case EXT_LZ4: {
			int n = LZ4_decompress_safe((const char*)in, (char*)out,
			                            (int)compLen, (int)uncompLen);
			if(n != (int)uncompLen) {
				std::cerr << "Error: LZ4_decompress_safe returned " << n
				          << ", expected " << uncompLen << std::endl;
				exit(1);
			}
			return;
		}
		case EXT_ZSTD: {
			size_t n = ZSTD_decompress(out, uncompLen, in, compLen);
			if(ZSTD_isError(n) || n != uncompLen) {
				std::cerr << "Error: ZSTD_decompress failed" << std::endl;
				exit(1);
			}
			return;
		}
#ifdef WITH_RANS
		case EXT_RANS0:
		case EXT_RANS1: {
			unsigned int outSize = (unsigned int)uncompLen;
			unsigned char* r = rans_uncompress_to_4x16(
				(unsigned char*)in, (unsigned int)compLen,
				(unsigned char*)out, &outSize);
			if(!r || outSize != uncompLen) {
				std::cerr << "Error: rans_uncompress_to_4x16 returned size "
				          << outSize << ", expected " << uncompLen << std::endl;
				exit(1);
			}
			return;
		}
#endif
		default:
			std::cerr << "Error: codec " << extCodecName(codec)
			          << " not available in this build" << std::endl;
			exit(1);
	}
}

// Worst-case compressed-output size estimate.
static size_t extCompressBound(ExtCodec codec, size_t inLen) {
	switch(codec) {
		case EXT_RAW:  return inLen;
		case EXT_LZ4:  return (size_t)LZ4_compressBound((int)inLen);
		case EXT_ZSTD: return ZSTD_compressBound(inLen);
#ifdef WITH_RANS
		case EXT_RANS0:
		case EXT_RANS1: {
			int order = (codec == EXT_RANS0) ? 0 : 1;
			return (size_t)rans_compress_bound_4x16((unsigned int)inLen, order);
		}
#endif
		default:       return inLen + 64 * 1024;
	}
}

// Sequential run-file writer. Buffers records into uncompressed blocks of
// ~kExtBlockTargetBytes, flushes one block at a time through the codec.
class ExtRunWriter {
public:
	ExtRunWriter(const char* path, ExtCodec codec,
	             const std::vector<std::string>& chrDict,
	             int numThreads = 1)
		: codec_(codec), numThreads_(numThreads),
		  totalRecords_(0), numBlocks_(0), numBlockRecords_(0)
	{
		fp_ = fopen(path, "wb");
		if(!fp_) { std::cerr << "Error: cannot create run file " << path
		                     << ": " << strerror(errno) << std::endl; exit(1); }
		uint8_t hdr[8];
		memcpy(hdr, "PSBR", 4);
		hdr[4] = 1;
		hdr[5] = (uint8_t)codec;
		hdr[6] = 0;
		hdr[7] = 0;
		fwrite(hdr, 1, 8, fp_);
		uint32_t chrCount = (uint32_t)chrDict.size();
		fwrite(&chrCount, sizeof(uint32_t), 1, fp_);
		for(const std::string& s : chrDict) {
			uint16_t len = (uint16_t)s.size();
			fwrite(&len, sizeof(uint16_t), 1, fp_);
			fwrite(s.data(), 1, len, fp_);
		}
		blockBuf_.reserve(kExtBlockTargetBytes + 1024);
	}

	void writeRecord(uint16_t chrIdx, int32_t beg, int32_t end,
	                 uint8_t strand, const char* tail, uint16_t tailLen)
	{
		const size_t need = 13 + (size_t)tailLen;
		if(!blockBuf_.empty() && blockBuf_.size() + need > kExtBlockTargetBytes)
			flushBlock();
		size_t off = blockBuf_.size();
		blockBuf_.resize(off + need);
		uint8_t* p = blockBuf_.data() + off;
		memcpy(p, &chrIdx, 2);
		memcpy(p + 2, &beg,  4);
		memcpy(p + 6, &end,  4);
		p[10] = strand;
		memcpy(p + 11, &tailLen, 2);
		if(tailLen) memcpy(p + 13, tail, tailLen);
		numBlockRecords_++;
		totalRecords_++;
	}

	void close() {
		flushBlock();
		uint8_t footer[20];
		memcpy(footer, "PSBE", 4);
		memcpy(footer + 4, &numBlocks_, sizeof(uint32_t));
		memcpy(footer + 8, &totalRecords_, sizeof(uint64_t));
		memset(footer + 16, 0, 4);
		fwrite(footer, 1, 20, fp_);
		fclose(fp_);
		fp_ = NULL;
	}

	~ExtRunWriter() { if(fp_) fclose(fp_); }

private:
	void flushBlock() {
		if(blockBuf_.empty()) return;
		uint32_t uncompLen = (uint32_t)blockBuf_.size();
		size_t maxComp = extCompressBound(codec_, uncompLen);
		std::vector<uint8_t> compBuf(maxComp);
		uint32_t compLen = (uint32_t)extCompressBlock(codec_,
			blockBuf_.data(), uncompLen, compBuf.data(), maxComp, numThreads_);
		uint32_t hdr[3] = { uncompLen, compLen, numBlockRecords_ };
		fwrite(hdr, sizeof(uint32_t), 3, fp_);
		fwrite(compBuf.data(), 1, compLen, fp_);
		blockBuf_.clear();
		numBlockRecords_ = 0;
		numBlocks_++;
	}

	FILE* fp_;
	ExtCodec codec_;
	int numThreads_;
	std::vector<uint8_t> blockBuf_;
	uint64_t totalRecords_;
	uint32_t numBlocks_;
	uint32_t numBlockRecords_;
};

// One-record-at-a-time reader. peek() returns the next record's key/data;
// advance() consumes it and lazily refills the block buffer when exhausted.
class ExtRunReader {
public:
	bool eof;
	uint16_t chrIdx;
	int32_t  beg;
	int32_t  end;
	uint8_t  strand;
	uint16_t tailLen;
	const uint8_t* tailPtr;
	std::vector<std::string> chrDict;

	ExtRunReader(const std::string& path)
		: eof(false), chrIdx(0), beg(0), end(0), strand('+'), tailLen(0), tailPtr(NULL),
		  blockPos_(0), blockEnd_(0), blocksRemaining_(0)
	{
		fp_ = fopen(path.c_str(), "rb");
		if(!fp_) { std::cerr << "Error: cannot open run file " << path
		                     << ": " << strerror(errno) << std::endl; exit(1); }
		uint8_t hdr[8];
		if(fread(hdr, 1, 8, fp_) != 8 || memcmp(hdr, "PSBR", 4) != 0) {
			std::cerr << "Error: bad run-file header in " << path << std::endl;
			exit(1);
		}
		if(hdr[4] != 1) {
			std::cerr << "Error: unsupported run-file version " << (int)hdr[4]
			          << " in " << path << std::endl;
			exit(1);
		}
		codec_ = (ExtCodec)hdr[5];
		uint32_t chrCount;
		fread(&chrCount, sizeof(uint32_t), 1, fp_);
		chrDict.reserve(chrCount);
		for(uint32_t i = 0; i < chrCount; i++) {
			uint16_t len;
			fread(&len, sizeof(uint16_t), 1, fp_);
			std::string s(len, '\0');
			fread(s.data(), 1, len, fp_);
			chrDict.push_back(std::move(s));
		}
		// Read footer to learn block count, then rewind to the first block.
		long blocksStart = ftell(fp_);
		fseek(fp_, -20, SEEK_END);
		uint8_t footer[20];
		fread(footer, 1, 20, fp_);
		if(memcmp(footer, "PSBE", 4) != 0) {
			std::cerr << "Error: bad run-file footer in " << path << std::endl;
			exit(1);
		}
		memcpy(&blocksRemaining_, footer + 4, sizeof(uint32_t));
		fseek(fp_, blocksStart, SEEK_SET);
		advance();
	}

	~ExtRunReader() { if(fp_) fclose(fp_); }

	void advance() {
		if(blockPos_ >= blockEnd_) {
			if(!loadNextBlock()) { eof = true; tailPtr = NULL; return; }
		}
		const uint8_t* p = blockBuf_.data() + blockPos_;
		memcpy(&chrIdx, p, 2);
		memcpy(&beg,  p + 2, 4);
		memcpy(&end,  p + 6, 4);
		strand = p[10];
		memcpy(&tailLen, p + 11, 2);
		tailPtr = p + 13;
		blockPos_ += 13 + tailLen;
	}

private:
	bool loadNextBlock() {
		if(blocksRemaining_ == 0) return false;
		uint32_t hdr[3];
		if(fread(hdr, sizeof(uint32_t), 3, fp_) != 3) return false;
		uint32_t uncompLen = hdr[0];
		uint32_t compLen   = hdr[1];
		blockBuf_.resize(uncompLen);
		if(uncompLen == compLen && codec_ == EXT_RAW) {
			fread(blockBuf_.data(), 1, uncompLen, fp_);
		} else {
			std::vector<uint8_t> compBuf(compLen);
			fread(compBuf.data(), 1, compLen, fp_);
			extDecompressBlock(codec_, compBuf.data(), compLen,
			                   blockBuf_.data(), uncompLen);
		}
		blockPos_ = 0;
		blockEnd_ = uncompLen;
		blocksRemaining_--;
		return true;
	}

	FILE* fp_;
	ExtCodec codec_;
	std::vector<uint8_t> blockBuf_;
	size_t blockPos_;
	size_t blockEnd_;
	uint32_t blocksRemaining_;
};

// ============================================================================
// External merge driver
// ============================================================================
//
// Pass 1: stream the input via fread; for each line, parse with
// parseBedLineFull; assign chrIdx (per-run dictionary, intern via map);
// store a per-line ExtEntry plus the tail bytes in a per-run Arena. When the
// arena+entry footprint hits the budget, sort entries with radixSort64 over
// packed (chrIdx<<32 | beg) keys and flush a run file. Repeat.
//
// Merge: open all run files; build a min-heap of (key, runIdx); pop, emit
// to stdout (reconstruct chr + write tab-separated text), advance that run.
//
// chrIdx assignment is per-run during pass 1. Each run ships its own chr
// dictionary in its header. The merge phase remaps each run's chr-name back
// to the global sorted dictionary built across all runs.
//
// The chr-name compare in the heap is done via a precomputed
// runChrToGlobal[runIdx][localChrIdx] -> globalChrIdx table, and the heap
// key is (globalChrIdx, beg) packed. Output emits each line as
// global_chr_name + "\t" + beg + "\t" + end + tail.

static int extMergeSort(const std::string& inputFile,
                        size_t memBudget,
                        ExtCodec codec,
                        const std::string& tmpDirOpt,
                        bool naturalSort,
                        int numThreads,
                        bool verbose)
{
	if(memBudget == 0) memBudget = 1ULL << 30;  // default 1 GB
	std::string tmpDir = tmpDirOpt;
	if(tmpDir.empty()) {
		const char* envTmp = getenv("TMPDIR");
		tmpDir = envTmp ? envTmp : "/tmp";
	}

	time_t tstart, tend;
	if(verbose) time(&tstart);

	FILE* in = fopen(inputFile.c_str(), "r");
	if(!in) {
		std::cerr << "Error: cannot open " << inputFile << ": "
		          << strerror(errno) << std::endl;
		return 1;
	}

	// Per-run state (resets each time we flush a run).
	std::vector<ExtEntry> entries;
	Arena* tailArena = new Arena(1UL << 20, 1UL << 24);
	std::unordered_map<std::string, uint16_t> chrIdxMap;
	std::vector<std::string> chrDict;
	size_t bytesUsed = 0;  // approximate accounting for entries + tails
	std::vector<std::string> runPaths;

	auto flushRun = [&]() {
		if(entries.empty()) return;
		size_t n = entries.size();

		// chrDict is in encounter order. Sort it (alphabetically or
		// naturally) and remap each entry's chrIdx so that the run's
		// local chrIdx order matches the global merge's chrIdx order.
		// Without this, a run that saw "chr10" before "chr1" would emit
		// chr10 records before chr1 records, breaking the k-way merge
		// invariant ("smallest local key == smallest global key").
		std::vector<uint32_t> sortedOrder(chrDict.size());
		for(uint32_t i = 0; i < chrDict.size(); i++) sortedOrder[i] = i;
		if(naturalSort) {
			std::sort(sortedOrder.begin(), sortedOrder.end(),
			          [&](uint32_t a, uint32_t b) {
			              return naturalChrLess(chrDict[a], chrDict[b]);
			          });
		} else {
			std::sort(sortedOrder.begin(), sortedOrder.end(),
			          [&](uint32_t a, uint32_t b) {
			              return chrDict[a] < chrDict[b];
			          });
		}
		std::vector<uint16_t> oldToNew(chrDict.size());
		std::vector<std::string> sortedDict(chrDict.size());
		for(uint32_t newIdx = 0; newIdx < sortedOrder.size(); newIdx++) {
			oldToNew[sortedOrder[newIdx]] = (uint16_t)newIdx;
			sortedDict[newIdx] = std::move(chrDict[sortedOrder[newIdx]]);
		}
		for(ExtEntry& e : entries) e.chrIdx = oldToNew[e.chrIdx];

		std::vector<uint64_t> keys(n);
		std::vector<uint32_t> order(n);
		for(size_t i = 0; i < n; i++) {
			keys[i]  = ((uint64_t)entries[i].chrIdx << 32) | (uint32_t)entries[i].beg;
			order[i] = (uint32_t)i;
		}
		radixSort64(keys.data(), order.data(), n, numThreads);

		char path[512];
		snprintf(path, sizeof(path), "%s/piosort.run.%u.%d.tmp",
		         tmpDir.c_str(), (unsigned)runPaths.size(), (int)getpid());
		ExtRunWriter writer(path, codec, sortedDict, numThreads);
		for(size_t i = 0; i < n; i++) {
			const ExtEntry& e = entries[order[i]];
			writer.writeRecord(e.chrIdx, e.beg, e.end, e.strand, e.tailPtr, e.tailLen);
		}
		writer.close();
		runPaths.emplace_back(path);

		entries.clear();
		entries.shrink_to_fit();
		delete tailArena;
		tailArena = new Arena(1UL << 20, 1UL << 24);
		chrIdxMap.clear();
		chrDict.clear();
		bytesUsed = 0;
	};

	// Read input line by line via getline (handles arbitrary line lengths).
	char* lineBuf = NULL;
	size_t lineCap = 0;
	ssize_t llen;
	bool sawHeader = false;
	while((llen = getline(&lineBuf, &lineCap, in)) != -1) {
		if(llen == 0) continue;
		if(lineBuf[llen - 1] == '\n') { lineBuf[--llen] = '\0'; }
		if(llen > 0 && lineBuf[llen - 1] == '\r') { lineBuf[--llen] = '\0'; }
		// Header passthrough (only at start of file, by convention).
		if(!sawHeader && (lineBuf[0] == '#' ||
		   strncmp(lineBuf, "track ",   6) == 0 ||
		   strncmp(lineBuf, "browser ", 8) == 0)) {
			fputs_unlocked(lineBuf, stdout);
			fputc_unlocked('\n', stdout);
			continue;
		}
		sawHeader = true;
		if(llen == 0) continue;

		const char* chrPtr;
		int chrLen;
		int beg = 0, end = 0;
		char weight[kWeightBufSize]; weight[0] = '0'; weight[1] = '\0';
		char strandChar = '+';
		const char* tailPtr = "";
		int n = parseBedLine3(lineBuf, &chrPtr, &chrLen, &beg, &end, &tailPtr);
		if(n < 3) {
			std::cerr << "Error in parsing line: " << lineBuf << std::endl;
			free(lineBuf); fclose(in); return 1;
		}
		(void)weight; (void)strandChar;

		// Intern chr name (per-run dictionary).
		std::string chrName(chrPtr, chrLen);
		uint16_t chrIdx;
		auto it = chrIdxMap.find(chrName);
		if(it != chrIdxMap.end()) {
			chrIdx = it->second;
		} else {
			if(chrDict.size() >= 65535) {
				std::cerr << "Error: more than 65535 distinct chromosomes per run"
				          << std::endl;
				return 1;
			}
			chrIdx = (uint16_t)chrDict.size();
			chrDict.push_back(chrName);
			chrIdxMap[chrName] = chrIdx;
		}

		// Append tail to the arena.
		size_t tailLen = strlen(tailPtr);
		if(tailLen > 65535) {
			std::cerr << "Error: BED line tail exceeds 64 KB" << std::endl;
			return 1;
		}
		const char* tailCopy = (tailLen > 0) ? tailArena->alloc(tailPtr, tailLen) : "";

		ExtEntry e;
		e.chrIdx  = chrIdx;
		e.beg     = beg;
		e.end     = end;
		e.strand  = (uint8_t)strandChar;
		e.tailLen = (uint16_t)tailLen;
		e.tailPtr = tailCopy;
		entries.push_back(e);

		bytesUsed += sizeof(ExtEntry) + tailLen + 2;
		if(bytesUsed >= memBudget) flushRun();
	}
	free(lineBuf);
	fclose(in);
	flushRun();
	delete tailArena;
	tailArena = NULL;

	if(verbose) {
		time(&tend);
		uint64_t totalTempBytes = 0;
		for(const std::string& p : runPaths) {
			struct stat st;
			if(stat(p.c_str(), &st) == 0) totalTempBytes += (uint64_t)st.st_size;
		}
		std::cerr << "Pass 1: " << runPaths.size() << " runs, "
		          << "temp_bytes=" << totalTempBytes
		          << " (" << ((double)totalTempBytes / (1024.0 * 1024.0)) << " MiB), "
		          << "codec=" << extCodecName(codec) << ", "
		          << (long)(tend - tstart) << " s" << std::endl;
		time(&tstart);
	}

	// ----- Merge phase -----
	// Open all runs; build global chr dictionary from union of per-run dicts;
	// remap each run's local chrIdx to global.
	std::vector<std::unique_ptr<ExtRunReader>> readers;
	readers.reserve(runPaths.size());
	for(const std::string& p : runPaths) {
		readers.emplace_back(std::make_unique<ExtRunReader>(p));
	}

	// Build global chr dictionary (sorted lex or natural).
	std::vector<std::string> globalChroms;
	{
		std::unordered_set<std::string> seen;
		for(auto& r : readers) {
			for(const std::string& c : r->chrDict) {
				if(seen.insert(c).second) globalChroms.push_back(c);
			}
		}
		if(naturalSort) std::sort(globalChroms.begin(), globalChroms.end(), naturalChrLess);
		else            std::sort(globalChroms.begin(), globalChroms.end());
	}
	std::unordered_map<std::string, uint32_t> globalIdxMap;
	for(uint32_t i = 0; i < globalChroms.size(); i++) globalIdxMap[globalChroms[i]] = i;

	// Per-run remap tables: run-local chrIdx -> global chrIdx.
	std::vector<std::vector<uint32_t>> remap(readers.size());
	for(size_t r = 0; r < readers.size(); r++) {
		remap[r].resize(readers[r]->chrDict.size());
		for(size_t i = 0; i < readers[r]->chrDict.size(); i++)
			remap[r][i] = globalIdxMap[readers[r]->chrDict[i]];
	}

	// Min-heap of (globalKey, runIdx). globalKey = (globalChrIdx<<32) | beg.
	struct HeapElem { uint64_t key; uint32_t runIdx; };
	auto heapCmp = [](const HeapElem& a, const HeapElem& b) {
		return a.key > b.key;  // min-heap
	};
	std::priority_queue<HeapElem, std::vector<HeapElem>, decltype(heapCmp)> heap(heapCmp);
	for(uint32_t r = 0; r < readers.size(); r++) {
		if(!readers[r]->eof) {
			uint64_t k = ((uint64_t)remap[r][readers[r]->chrIdx] << 32)
			           | (uint32_t)readers[r]->beg;
			heap.push({k, r});
		}
	}

	// Output buffer (8 MiB) is set up by main(); we reuse stdout.
	char numBuf[32];
	while(!heap.empty()) {
		HeapElem top = heap.top(); heap.pop();
		ExtRunReader* rd = readers[top.runIdx].get();
		const std::string& chr = globalChroms[(uint32_t)(top.key >> 32)];
		fputs_unlocked(chr.c_str(), stdout);
		fputc_unlocked('\t', stdout);
		int n = writeUInt(numBuf, rd->beg);
		numBuf[n++] = '\t';
		n += writeUInt(numBuf + n, rd->end);
		fwrite_unlocked(numBuf, 1, n, stdout);
		if(rd->tailLen) fwrite_unlocked(rd->tailPtr, 1, rd->tailLen, stdout);
		fputc_unlocked('\n', stdout);
		rd->advance();
		if(!rd->eof) {
			uint64_t k = ((uint64_t)remap[top.runIdx][rd->chrIdx] << 32)
			           | (uint32_t)rd->beg;
			heap.push({k, top.runIdx});
		}
	}

	// Cleanup.
	readers.clear();
	for(const std::string& p : runPaths) unlink(p.c_str());

	if(verbose) {
		time(&tend);
		std::cerr << "Merge: " << (long)(tend - tstart) << " s" << std::endl;
	}
	return 0;
}

// ============================================================================
// MULTI-PASS NO-WRITES PATH (--multi-pass)
// ----------------------------------------------------------------------------
// In-RAM sort for inputs ~1-5x RAM, with zero disk writes (SSD-wear-friendly
// fallback). Pass 1 streams the file via getline and builds a histogram
// keyed by (chrIdx, beg >> kPromBucketShift), tracking byte size per bucket.
// The histogram is then bin-packed into K consecutive groups, each <= the
// memory budget. Passes 2..K+1 re-stream the input; each pass keeps only
// records that fall in its group's bucket range, sorts them with
// radixSort64, and emits to stdout.
//
// Total disk reads = (K+1) x file_size, total disk writes = 0. K = ceil(
// file_size / budget). For "medium" inputs (1-5x RAM) this is 2-6 reads,
// dwarfed by what bedops external-merge would write to NAND. Above ~5x RAM
// the K-pass cost becomes quadratic; --external-merge is asymptotically
// better there.
//
// Pathological inputs (a single 1 MB position bucket bigger than the
// budget) error out with a pointer at --max-mem / --external-merge.
// ============================================================================

constexpr int kPromBucketShift = 20;  // 1 MB position quantum per bucket

struct PromHistEntry {
	uint64_t bytes;
	uint64_t records;
};

static int multiPassSort(const std::string& inputFile,
                         size_t budget,
                         bool naturalSort,
                         int numThreads,
                         bool verbose)
{
	if(budget == 0) budget = 1ULL << 30;

	time_t tstart, tend;
	if(verbose) time(&tstart);

	// ----- Pass 1: histogram -----
	FILE* in = fopen(inputFile.c_str(), "r");
	if(!in) {
		std::cerr << "Error: cannot open " << inputFile << ": "
		          << strerror(errno) << std::endl;
		return 1;
	}

	std::unordered_map<std::string, uint32_t> chrEncMap;
	std::vector<std::string> chrNamesEnc;
	std::unordered_map<uint64_t, PromHistEntry> hist;

	char* lineBuf = NULL;
	size_t lineCap = 0;
	ssize_t llen;
	bool sawHeader = false;
	uint64_t totalLines = 0;
	uint64_t totalBytes = 0;

	while((llen = getline(&lineBuf, &lineCap, in)) != -1) {
		if(llen == 0) continue;
		size_t sz = (size_t)llen;
		if(lineBuf[sz-1] == '\n') { lineBuf[--sz] = '\0'; }
		if(sz > 0 && lineBuf[sz-1] == '\r') { lineBuf[--sz] = '\0'; }
		if(sz == 0) continue;
		if(!sawHeader && (lineBuf[0] == '#' ||
		    strncmp(lineBuf, "track ",   6) == 0 ||
		    strncmp(lineBuf, "browser ", 8) == 0)) {
			fputs_unlocked(lineBuf, stdout);
			fputc_unlocked('\n', stdout);
			continue;
		}
		sawHeader = true;

		const char* chrPtr;
		int chrLen;
		int beg = 0, end = 0;
		const char* tailPtr = "";
		if(parseBedLine3(lineBuf, &chrPtr, &chrLen, &beg, &end, &tailPtr) < 3) {
			std::cerr << "Error in parsing line: " << lineBuf << std::endl;
			free(lineBuf); fclose(in); return 1;
		}
		std::string chrName(chrPtr, chrLen);
		uint32_t encIdx;
		auto it = chrEncMap.find(chrName);
		if(it != chrEncMap.end()) {
			encIdx = it->second;
		} else {
			encIdx = (uint32_t)chrNamesEnc.size();
			chrNamesEnc.push_back(chrName);
			chrEncMap[chrName] = encIdx;
		}
		uint64_t bid = ((uint64_t)encIdx << 32) | (uint32_t)(beg >> kPromBucketShift);
		PromHistEntry& e = hist[bid];
		e.bytes   += (uint64_t)sz + 1;  // include line terminator
		e.records += 1;
		totalLines += 1;
		totalBytes += (uint64_t)sz + 1;
	}
	free(lineBuf);
	fclose(in);

	// Sort chr names alphabetically (or naturally) to match the global emit order.
	std::vector<uint32_t> sortOrder(chrNamesEnc.size());
	for(uint32_t i = 0; i < sortOrder.size(); i++) sortOrder[i] = i;
	if(naturalSort) std::sort(sortOrder.begin(), sortOrder.end(),
	    [&](uint32_t a, uint32_t b) { return naturalChrLess(chrNamesEnc[a], chrNamesEnc[b]); });
	else std::sort(sortOrder.begin(), sortOrder.end(),
	    [&](uint32_t a, uint32_t b) { return chrNamesEnc[a] < chrNamesEnc[b]; });

	std::vector<uint32_t> encToSorted(chrNamesEnc.size());
	for(uint32_t i = 0; i < sortOrder.size(); i++)
		encToSorted[sortOrder[i]] = i;

	// Remap histogram keys (encChr -> sortedChr) and sort buckets in
	// global sort order.
	struct Bucket { uint64_t bid; uint64_t bytes; uint64_t records; };
	std::vector<Bucket> buckets;
	buckets.reserve(hist.size());
	for(auto& kv : hist) {
		uint32_t encIdx = (uint32_t)(kv.first >> 32);
		uint64_t newBid = ((uint64_t)encToSorted[encIdx] << 32) | (uint32_t)kv.first;
		buckets.push_back({newBid, kv.second.bytes, kv.second.records});
	}
	std::sort(buckets.begin(), buckets.end(),
	          [](const Bucket& a, const Bucket& b) { return a.bid < b.bid; });

	// ----- Bin-pack into K groups -----
	struct Group { uint64_t startBid; uint64_t endBid; uint64_t bytes; uint64_t records; };
	std::vector<Group> groups;
	if(!buckets.empty()) {
		// Per-record overhead estimate for ExtEntry + arena fragmentation:
		// ~20 bytes / record beyond the line itself. Account for this so
		// the group's true peak RAM stays within budget.
		const uint64_t perRecOverhead = 20;
		uint64_t curStart = buckets[0].bid;
		uint64_t curBytes = 0;
		uint64_t curRecords = 0;
		for(size_t i = 0; i < buckets.size(); i++) {
			uint64_t bytesWithOverhead = buckets[i].bytes
			                           + buckets[i].records * perRecOverhead;
			if(bytesWithOverhead > budget) {
				std::cerr << "Error: a single 1 MB-quantum bucket holds "
				          << buckets[i].records << " records / "
				          << buckets[i].bytes << " bytes (with overhead "
				          << bytesWithOverhead << "), exceeding --max-mem "
				          << budget << ". Raise --max-mem or use --external-merge."
				          << std::endl;
				return 1;
			}
			if(curBytes + bytesWithOverhead > budget) {
				groups.push_back({curStart, buckets[i].bid, curBytes, curRecords});
				curStart = buckets[i].bid;
				curBytes = 0;
				curRecords = 0;
			}
			curBytes += bytesWithOverhead;
			curRecords += buckets[i].records;
		}
		// Close last group (use last bucket's bid + 1 as exclusive end).
		uint64_t endBid = buckets.back().bid + 1;
		groups.push_back({curStart, endBid, curBytes, curRecords});
	}

	if(verbose) {
		time(&tend);
		std::cerr << "Pass 1: " << totalLines << " records, "
		          << buckets.size() << " buckets, "
		          << groups.size() << " groups (K), "
		          << "total_bytes=" << totalBytes
		          << " ("
		          << ((double)totalBytes / (1024.0 * 1024.0 * 1024.0))
		          << " GiB), "
		          << (long)(tend - tstart) << " s" << std::endl;
		time(&tstart);
	}

	// ----- Passes 2..K+1: re-stream input per group, sort, emit -----
	std::vector<ExtEntry> entries;
	Arena* arena = new Arena(1UL << 20, 1UL << 24);

	for(size_t gi = 0; gi < groups.size(); gi++) {
		Group& g = groups[gi];
		entries.clear();
		entries.reserve(g.records);
		delete arena;
		arena = new Arena(1UL << 20, 1UL << 24);

		in = fopen(inputFile.c_str(), "r");
		if(!in) {
			std::cerr << "Error: cannot reopen " << inputFile << std::endl;
			delete arena; return 1;
		}

		lineBuf = NULL; lineCap = 0;
		bool seenBody = false;
		while((llen = getline(&lineBuf, &lineCap, in)) != -1) {
			if(llen == 0) continue;
			size_t sz = (size_t)llen;
			if(lineBuf[sz-1] == '\n') { lineBuf[--sz] = '\0'; }
			if(sz > 0 && lineBuf[sz-1] == '\r') { lineBuf[--sz] = '\0'; }
			if(sz == 0) continue;
			if(!seenBody && (lineBuf[0] == '#' ||
			    strncmp(lineBuf, "track ",   6) == 0 ||
			    strncmp(lineBuf, "browser ", 8) == 0)) {
				continue;  // already emitted in pass 1
			}
			seenBody = true;

			const char* chrPtr;
			int chrLen;
			int beg = 0, end = 0;
			const char* tailPtr = "";
			if(parseBedLine3(lineBuf, &chrPtr, &chrLen, &beg, &end, &tailPtr) < 3) {
				std::cerr << "Error in parsing line: " << lineBuf << std::endl;
				free(lineBuf); fclose(in); delete arena; return 1;
			}
			std::string chrName(chrPtr, chrLen);
			auto it = chrEncMap.find(chrName);
			if(it == chrEncMap.end()) {
				std::cerr << "Error: chromosome '" << chrName
				          << "' appeared in pass " << (gi + 2)
				          << " but not in pass 1 (input changed during sort?)"
				          << std::endl;
				free(lineBuf); fclose(in); delete arena; return 1;
			}
			uint32_t sortedIdx = encToSorted[it->second];
			uint64_t bid = ((uint64_t)sortedIdx << 32)
			             | (uint32_t)(beg >> kPromBucketShift);
			if(bid < g.startBid || bid >= g.endBid) continue;

			// Match: copy full line text into arena.
			const char* lineCopy = arena->alloc(lineBuf, sz);
			ExtEntry e;
			e.chrIdx  = (uint16_t)sortedIdx;
			e.beg     = beg;
			e.end     = end;
			e.strand  = '+';
			e.tailLen = (uint16_t)sz;     // re-use tailLen as "full line length"
			e.tailPtr = lineCopy;          // tailPtr points at the entire line
			entries.push_back(e);
		}
		free(lineBuf);
		fclose(in);

		// Sort group entries by (chrIdx, beg) and emit.
		size_t n = entries.size();
		std::vector<uint64_t> keys(n);
		std::vector<uint32_t> order(n);
		for(size_t i = 0; i < n; i++) {
			keys[i]  = ((uint64_t)entries[i].chrIdx << 32) | (uint32_t)entries[i].beg;
			order[i] = (uint32_t)i;
		}
		radixSort64(keys.data(), order.data(), n, numThreads);
		for(size_t i = 0; i < n; i++) {
			const ExtEntry& e = entries[order[i]];
			fwrite_unlocked(e.tailPtr, 1, e.tailLen, stdout);
			fputc_unlocked('\n', stdout);
		}

		if(verbose) {
			time(&tend);
			std::cerr << "Pass " << (gi + 2) << "/" << (groups.size() + 1)
			          << " (group " << (gi + 1) << "/" << groups.size()
			          << "): " << n << " records emitted in "
			          << (long)(tend - tstart) << " s" << std::endl;
			time(&tstart);
		}
	}

	delete arena;
	return 0;
}

// ============================================================================
// CLI / MAIN
// ----------------------------------------------------------------------------
// CLI11 option parsing, input dispatch (file mmap / stdin slurp / gzip
// slurp), and routing to lowMemSortMmap or the classic + bucket sort path.
// Owns parseMemSize (for --max-mem) and the slurpStream helper.
// ============================================================================

// Parse a memory-size string like "4G", "500M", "1024K", or a bare byte count.
// Returns 0 on empty/malformed input. Suffixes use binary multipliers (1024).
static size_t parseMemSize(const std::string& s)
{
	if(s.empty()) return 0;
	size_t mult = 1;
	std::string num = s;
	char back = s.back();
	if(back == 'G' || back == 'g') { mult = 1ULL << 30; num.pop_back(); }
	else if(back == 'M' || back == 'm') { mult = 1ULL << 20; num.pop_back(); }
	else if(back == 'K' || back == 'k') { mult = 1ULL << 10; num.pop_back(); }
	if(num.empty()) return 0;
	try
	{
		long long v = std::stoll(num);
		if(v < 0) return 0;
		return (size_t)v * mult;
	}
	catch(...) { return 0; }
}

int main(int argc, char *argv[])
{
	// 8MB output buffer — amortises write syscalls on 100M+ read files
	char* outBuf = (char*) malloc(1 << 23);
	setvbuf(stdout, outBuf, _IOFBF, 1 << 23);

	CLI::App app{"Ultra fast bed file sorter\nPiotr Balwierz, 2012-2026\n\n"
		"Results are equivalent to \"LC_ALL=C sort -k1,1 -k2,2n file.bed\"\n"
		"or \"sort -k1,1 -k2,2n -k3,3n file.bed\" if --sort=b enabled.\n"
		"Defaults to all cores; pass -t 1 for single-threaded.\n"
		"For files > ~1 M reads, --low-mem-ssd is the recommended fast path.\n\n"
		"Input file should contain a new line character in the end of the last line.\n"
		"Gzip-compressed input (.gz) is transparently decompressed via gzip.\n"
		"BED header lines (track/browser/#) are passed through unchanged.\n\n"
		"Compilation-time limits:\n"
		"  Line length: no fixed limit (stdin/gzip uses getline; mmap uses memchr).\n"
		"  Chromosome name: no fixed limit (stored as pointer+length).\n"
		"  Weight field (BED col 4, --collapse): " + to_string(kWeightBufSize - 1) + " B; over-long values are rejected with an error.\n"
		"  Read number limit: " + to_string((unsigned long long)UINT32_MAX - 1) + " (all paths).\n"};
	app.set_version_flag("-V,--version", VERSION_STRING);

	string inputFile;
	char sortMode = 's';
	int fCollapse = 0;
	int lowMemSSD = 0;
	int extMerge = 0;
	int multiPass = 0;
	// zstd default: fastest compressed codec across all builds, best
	// compression ratio (~0.33x of text input on a synthetic BED6, ~0.45x
	// of binary-uncompressed temps), and at 200 M records / 256 M budget
	// it is *faster* than `raw` end-to-end (compression saves more disk
	// I/O than its CPU costs). See benchmark/bench_extmerge.sh.
	std::string extCodecStr = "zstd";
	std::string extTmpDir;
	int numThreads = 0;
	bool naturalSort = false;
	std::string maxMemStr;
	std::string outputFile;
	bool verbose = false;

	app.add_option("input-file", inputFile,
		"input file; \"-\" reads from stdin; .gz files are decompressed automatically")
		->required();
	app.add_option("-o,--output", outputFile,
		"write sorted output to this file (default: stdout)");
	app.add_option("-s,--sort", sortMode,
		"s: start coordinate [default]\n"
		"b: sort also by end coordinate\n"
		"5: sort by 5' end (respects strand)")
		->default_val('s');
	app.add_flag("--collapse", fCollapse,
		"collapse BEDWEIGHT by summing weights of the reads and truncate the coordinates, "
		"set strand to \"+\", id to \".\". Makes no sense for --sort=b");
	app.add_flag("--low-mem-ssd", lowMemSSD,
		"low-memory two-pass file mode (SSD-friendly): keeps line offsets in RAM and "
		"sorts one chromosome at a time; slower than default but uses less peak RAM");
	app.add_flag("--external-merge", extMerge,
		"external merge sort: stream the input, build sorted runs in RAM up to "
		"--max-mem (default 1G), spill each run as a binary temp file, k-way "
		"merge to stdout. Use for inputs > RAM. Temp files go to --tmpdir or "
		"$TMPDIR. (only the default coordinate sort is supported)");
	app.add_flag("--multi-pass", multiPass,
		"K-pass scan, no disk writes: pass 1 builds a histogram + bin-packs "
		"into K consecutive (chr,beg) groups <= --max-mem; passes 2..K+1 "
		"re-stream the input filtering by group, sort + emit. Use for "
		"medium inputs (~1-5x RAM) where SSD-wear matters and you don't "
		"want temp files. K = ceil(file_size / --max-mem). "
		"(only the default coordinate sort; file input only)");
	app.add_option("--merge-codec", extCodecStr,
		"compression codec for --external-merge temp files: raw|lz4|zstd"
		" (rans0|rans1 also available in WITH_BAM builds). zstd is the "
		"default — fastest compressed codec, best ratio, and at large "
		"sizes it is faster than `raw` end-to-end because compression "
		"saves more disk I/O than its CPU costs.")
		->default_val("zstd");
	app.add_option("--tmpdir", extTmpDir,
		"temp directory for --external-merge run files (default: $TMPDIR or /tmp)");
	app.add_option("-t,--threads", numThreads,
		"number of threads for the classic sort path "
		"(0 = all cores; 1 = single-threaded std::sort)")
		->default_val(0);
	app.add_flag("-n,--natural-sort", naturalSort,
		"sort chromosomes in natural order (chr2 < chr10) instead of "
		"lexicographic order (chr10 < chr2)");
	app.add_option("--max-mem", maxMemStr,
		"memory cap on per-chromosome scratch in --low-mem-ssd parallel "
		"emit (e.g. 4G, 500M, 2048K, or bare bytes). Without this, "
		"--low-mem-ssd runs uncapped. Set only to PREVENT OOM, not to "
		"optimise memory.");
	app.add_flag("-v,--verbose", verbose,
		"print parsing / sorting timing and per-chromosome length info to stderr "
		"(default: silent)");

	CLI11_PARSE(app, argc, argv);
	const size_t maxMemBytes = parseMemSize(maxMemStr);

	if(numThreads <= 0)
	{
		unsigned hw = std::thread::hardware_concurrency();
		numThreads = hw ? (int)hw : 1;
	}
	// TBB owns the worker pool used by std::execution::par. Capping it here
	// gives `--threads N` the same behavior the OpenMP build had. Set early
	// so it also applies to the BAM path's radixSort64 invocation.
	tbb::global_control tbbThreadCap(
		tbb::global_control::max_allowed_parallelism, (size_t)numThreads);

	// BAM input dispatch (opt-in at build time via WITH_BAM). Detect by .bam
	// extension; route to the in-RAM bamSortAndEmit path which reuses
	// radixSort64 and writes a coordinate-sorted BAM via htslib. None of the
	// BED-specific machinery below (stdout freopen, mmap, lowMemSortMmap,
	// bucket sort, etc.) applies — the BAM path opens outputFile directly via
	// sam_open instead of going through stdout.
#ifdef WITH_BAM
	if(inputFile.size() > 4 &&
	   (inputFile.compare(inputFile.size() - 4, 4, ".bam") == 0 ||
	    inputFile.compare(inputFile.size() - 4, 4, ".BAM") == 0))
	{
		if(fCollapse || lowMemSSD || naturalSort || sortMode != 's')
		{
			cerr << "Error: BAM input only supports default coordinate sort; "
			        "--collapse / --low-mem-ssd / --natural-sort / "
			        "--sort=b|5 are not implemented for BAM" << endl;
			return 1;
		}
		return bamSortAndEmit(inputFile, outputFile, numThreads, verbose);
	}
#endif

	// Multi-pass no-writes dispatch. Streams the input via getline K times
	// (1 histogram pass + K group passes); never writes to disk. For inputs
	// in the ~1-5x RAM range where SSD-wear avoidance matters more than
	// the (K+1)x read cost.
	if(multiPass)
	{
		if(fCollapse || lowMemSSD || extMerge || sortMode != 's')
		{
			cerr << "Error: --multi-pass only supports the default coordinate "
			        "sort (--sort=s) and is mutually exclusive with --collapse, "
			        "--low-mem-ssd, and --external-merge" << endl;
			return 1;
		}
		if(inputFile == "-") {
			cerr << "Error: --multi-pass requires a seekable input file "
			        "(stdin can only be read once)." << endl;
			return 1;
		}
		if(!outputFile.empty()) {
			if(!freopen(outputFile.c_str(), "w", stdout)) {
				cerr << "Error: cannot open " << outputFile << " for writing" << endl;
				return 1;
			}
		}
		return multiPassSort(inputFile, maxMemBytes, naturalSort, numThreads, verbose);
	}

	// External merge sort dispatch. Fires before mmap/slurp and bypasses the
	// normal input pipeline — extMergeSort opens the input file directly via
	// fread and writes its own output stream.
	if(extMerge)
	{
		if(fCollapse || lowMemSSD || sortMode != 's')
		{
			cerr << "Error: --external-merge only supports the default "
			        "coordinate sort (--sort=s); --collapse / --low-mem-ssd / "
			        "--sort=b|5 are not yet implemented in this path" << endl;
			return 1;
		}
		if(inputFile == "-") {
			cerr << "Error: --external-merge does not yet support stdin "
			        "input (would defeat the streaming-from-disk design); "
			        "pass a file path." << endl;
			return 1;
		}
		ExtCodec codec;
		if(parseExtCodec(extCodecStr, &codec) != 0) {
			cerr << "Error: unknown --merge-codec value '" << extCodecStr
			     << "' (valid: raw|lz4|zstd|rans0|rans1)" << endl;
			return 1;
		}
		if(!outputFile.empty()) {
			if(!freopen(outputFile.c_str(), "w", stdout)) {
				cerr << "Error: cannot open " << outputFile << " for writing" << endl;
				return 1;
			}
		}
		return extMergeSort(inputFile, maxMemBytes, codec, extTmpDir,
		                    naturalSort, numThreads, verbose);
	}

	// -o/--output redirects stdout to the named file. All downstream data
	// emit (fputs_unlocked, fwrite_unlocked, printf in collapse mode) goes
	// through stdout, so freopen-ing it here is a single-line redirect with
	// no other code changes needed. Diagnostics stay on stderr.
	if(!outputFile.empty())
	{
		if(!freopen(outputFile.c_str(), "w", stdout))
		{
			cerr << "Error: cannot open " << outputFile << " for writing" << endl;
			return 1;
		}
	}

	size_t readCount = 1;	// count == 0 is used as an end of a list.
	size_t currMaxReads = 1024;
	long fileSize = -1;

	char* mmapBase = NULL;
	size_t mmapSize = 0;
	bool useMmap = false;
	bool isGzip = false;
	FILE* fh = NULL;
	Arena* arena = NULL;

	// Helper: slurp an entire FILE* into a malloc'd buffer + 1 NUL byte.
	// Buffer doubles on growth. Returns NULL on OOM (caller errors out).
	auto slurpStream = [](FILE* in, size_t* outFill) -> char* {
		size_t cap = 1UL << 24;  // 16 MB initial
		char* buf = (char*) malloc(cap + 1);
		if(!buf) return NULL;
		size_t fill = 0;
		size_t n;
		while((n = fread(buf + fill, 1, cap - fill, in)) > 0)
		{
			fill += n;
			if(fill == cap)
			{
				cap *= 2;
				char* nb = (char*) realloc(buf, cap + 1);
				if(!nb) { free(buf); return NULL; }
				buf = nb;
			}
		}
		buf[fill] = '\0';
		*outFill = fill;
		return buf;
	};

	if(inputFile == "-")
	{
		// Slurp stdin into a buffer and feed it to the mmap-style parser.
		// This gets us zero-copy line pointers (no per-line memcpy into an
		// arena) at the cost of holding the input in memory.
		size_t fill = 0;
		mmapBase = slurpStream(stdin, &fill);
		if(!mmapBase) { cerr << "Error: out of memory slurping stdin" << endl; return 1; }
		mmapSize = fill;
		useMmap = true;
		const size_t chunkSize = 64UL * 1024 * 1024;
		arena = new Arena(chunkSize, chunkSize);  // for --collapse weight strings
	}
	else
	{
		// Detect gzip input by .gz extension
		isGzip = inputFile.size() > 3 &&
		         inputFile.compare(inputFile.size() - 3, 3, ".gz") == 0;

		if(isGzip)
		{
			// Build shell command with single-quote-safe filename escaping
			std::string cmd = "gzip -dc -- '";
			for(char c : inputFile)
			{
				if(c == '\'') cmd += "'\\''";
				else cmd += c;
			}
			cmd += "'";
			fh = popen(cmd.c_str(), "r");
			if(!fh)
			{
				cerr << "Error: cannot decompress " << inputFile << endl;
				return 1;
			}
			// Slurp the decompressed output into a buffer (same trick as stdin).
			size_t fill = 0;
			mmapBase = slurpStream(fh, &fill);
			pclose(fh);
			fh = NULL;
			isGzip = false;  // no longer need pclose at end
			if(!mmapBase) { cerr << "Error: out of memory slurping gzip stream" << endl; return 1; }
			mmapSize = fill;
			useMmap = true;
			const size_t chunkSize = 64UL * 1024 * 1024;
			arena = new Arena(chunkSize, chunkSize);  // for --collapse weight strings
		}
		else
		{
			int fd = open(inputFile.c_str(), O_RDONLY);
			if(fd < 0)
			{
				cerr << "Error opening " << inputFile << endl;
				return 1;
			}
			struct stat st;
			fstat(fd, &st);
			fileSize = st.st_size;

			if(fileSize == 0)
			{
				close(fd);
				return 0;
			}

			// Map one extra byte so the NUL sentinel past EOF is always accessible,
			// even when file size is an exact multiple of the page size.
			//
			// MAP_POPULATE pre-faults all pages at mmap() time, avoiding millions
			// of page-fault traps during the parse pass on cold-cache files.
			// MADV_HUGEPAGE hints the kernel to back the region with transparent
			// huge pages, cutting TLB pressure on large genomic files.
			// Both flags are Linux extensions; gated on availability so the
			// source still compiles on non-Linux systems if anyone ports.
			int mmapFlags = MAP_PRIVATE;
#ifdef MAP_POPULATE
			mmapFlags |= MAP_POPULATE;
#endif
			mmapBase = (char*) mmap(NULL, (size_t)fileSize + 1, PROT_READ | PROT_WRITE, mmapFlags, fd, 0);
			if(mmapBase == MAP_FAILED)
			{
				cerr << "Error: mmap failed for " << inputFile << endl;
				close(fd);
				return 1;
			}
			close(fd);
			madvise(mmapBase, (size_t)fileSize, MADV_SEQUENTIAL);
#ifdef MADV_HUGEPAGE
			madvise(mmapBase, (size_t)fileSize, MADV_HUGEPAGE);
#endif
			mmapSize = (size_t)fileSize;
			useMmap = true;

			// Estimate read count from first 64KB
			size_t scanSize = mmapSize < 65536 ? mmapSize : 65536;
			size_t nLines = 0;
			for(size_t i = 0; i < scanSize; i++)
				if(mmapBase[i] == '\n') nLines++;
			if(nLines > 0)
				currMaxReads = (size_t)((double)mmapSize / (double)scanSize * (double)nLines * 1.1);
			if(currMaxReads < 1024)
				currMaxReads = 1024;

			if(fCollapse)
			{
				size_t arenaInit = (size_t)(fileSize / 5);
				if(arenaInit < 4096) arenaInit = 4096;
				arena = new Arena(arenaInit, 64UL * 1024 * 1024);
			}
		}
	}

	if(lowMemSSD)
	{
		// All input paths land here with useMmap=true (file: native mmap;
		// stdin/gzip: slurped into a buffer and presented as if mmap'd).
		return lowMemSortMmap(mmapBase, mmapSize, fCollapse, sortMode,
		                      numThreads, naturalSort, maxMemBytes, verbose);
	}

	time_t tstart, tend;
	if(verbose) time(&tstart);

	string2chrInfoT chrInfo;
	seqread* reads = NULL;

	// Per-chunk arenas for parallel mmap + --collapse. Owned by main; freed at exit.
	std::vector<Arena*> chunkArenas;

	if(useMmap)
	{
		// Parallel parser dispatches to parseMmapSerial when N==1, or to
		// parseChunkMmap-per-chunk otherwise (incl. --collapse, with one Arena
		// per chunk so weight-string copies don't serialise on a shared mutex).
		// Header lines (track/browser/#) are emitted to stdout in a leading
		// pre-pass before chunking.
		size_t parsed = 0;
		reads = parseMmapDispatch(mmapBase, mmapSize,
		                          fCollapse, sortMode,
		                          numThreads, arena, &chunkArenas,
		                          chrInfo, parsed);
		readCount = parsed + 1;  // matches the legacy "readCount-1 = total reads" semantics
	}
	else
	{
		reads = (seqread*) malloc(currMaxReads * sizeof(seqread));
		pair<string2chrInfoT::iterator,bool> insResult =
			chrInfo.insert(make_pair(string(""), chrInfoT()));
		string2chrInfoT::iterator thisChrIt = insResult.first;
		parseLines<false>(mmapBase, mmapSize, fh,
		                  reads, readCount, currMaxReads,
		                  chrInfo, thisChrIt, arena,
		                  fCollapse, sortMode);
		chrInfo.erase("");
	}

	if(isGzip && fh) pclose(fh);
	else if(fh && fh != stdin) fclose(fh);

	if(verbose)
	{
		time(&tend);
		fprintf(stderr, "Reading has taken %d seconds\n",
		        (int)(difftime(tend, tstart)+0.5));
	}

	/***********************************
	 *  the actual sorting starts here *
	 * *********************************/
	std::vector<std::string> chroms;
	for(string2chrInfoT::iterator it=chrInfo.begin(); it!=chrInfo.end(); it++)
	{
		if(verbose) cerr << it->first << endl;
		chroms.push_back(it->first);
	}

	if(naturalSort)
		std::sort(chroms.begin(), chroms.end(), naturalChrLess);
	else
		std::sort(chroms.begin(), chroms.end());

	// Assign sequential chromosome indices so the classic sort path can compare
	// chromosomes with a single int compare instead of a string compare.
	for(int ci = 0; ci < (int)chroms.size(); ci++)
		chrInfo.find(chroms[ci])->second.idx = ci;

	size_t totalReads = readCount - 1;
	if(verbose)
	{
		cerr << "We have " << totalReads << " regions, "
		     << chroms.size() << " chromosomes." << endl;
		time(&tstart);
	}

	/*************************************************************
	 * CLASSIC SORT PATH — index-array sort with packed-key radix
	 *
	 * Sort an index array (4-byte ints) rather than 24-byte seqread structs.
	 * Steps:
	 *  1. Walk linked lists, stamp chrIdx (reusing next's storage), build index array.
	 *  2. radixSort64 / std::sort the index array.
	 *  3. Linear scan to print (and optionally collapse).
	 *************************************************************/

	uint32_t* order = (uint32_t*) malloc(totalReads * sizeof(uint32_t));
	// chromStart[ci] = start offset in order[] for chromosome ci. Used by
	// the --sort b per-chrom radix path; harmless overhead otherwise.
	std::vector<size_t> chromStart(chroms.size() + 1, 0);
	size_t oi = 0;
	for(uint32_t ci = 0; ci < (uint32_t)chroms.size(); ci++)
	{
		chromStart[ci] = oi;
		string2chrInfoT::iterator cit = chrInfo.find(chroms[ci]);
		uint32_t curr = cit->second.lastRead;
		while(curr)
		{
			uint32_t nxt = reads[curr].next;
			reads[curr].chrIdx = ci;
			order[oi++] = curr;
			curr = nxt;
		}
	}
	chromStart[chroms.size()] = oi;

	// --sort s and --sort 5 use a global radix sort over packed
	// (chrIdx, pos) 64-bit keys. --sort b's three keys (chrIdx, beg, end)
	// don't fit in 64 bits without lossy squeezing, so it goes through
	// the per-chrom radix path which sorts each chromosome's slice of
	// order[] independently with a (beg, end) 64-bit key.
	if(sortMode == 'b')
		sortIndicesPerChromB(order, chromStart, reads, numThreads);
	else
		sortIndicesDispatch(order, totalReads, reads, sortMode, numThreads);

	if(fCollapse)
	{
		size_t i = 0;
		while(i < totalReads)
		{
			uint32_t ri = order[i];
			uint32_t ci = reads[ri].chrIdx;
			int pos = reads[ri].beg;
			float sum = 0.0f;
			while(i < totalReads)
			{
				uint32_t rj = order[i];
				if(reads[rj].chrIdx != ci || reads[rj].beg != pos) break;
				sum += parseWeight(reads, rj);
				i++;
			}
			printf("%s\t%d\t%d\t.\t%g\t+\n", chroms[ci].c_str(), pos, pos+1, sum);
		}
	}
	else if(useMmap)
	{
		for(size_t i = 0; i < totalReads; i++)
		{
			fputs_unlocked(reads[order[i]].line, stdout);
			fputc_unlocked('\n', stdout);
		}
	}
	else
	{
		for(size_t i = 0; i < totalReads; i++)
		{
			uint32_t ri = order[i];
			writeBedLine(chroms[reads[ri].chrIdx].c_str(),
			             reads[ri].beg, reads[ri].end,
			             reads[ri].line, stdout);
		}
	}
	free(order);

	if(verbose)
	{
		time(&tend);
		fprintf(stderr, "Sorting has taken %d seconds\n",
		        (int)(difftime(tend, tstart)+0.5));
	}
	return 0;
}
