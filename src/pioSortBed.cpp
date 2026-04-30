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
#include <mutex>
#include <thread>
#include <time.h>
#include <tbb/global_control.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "CLI11.hpp"

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
// chrInfoT:   per-chromosome bookkeeping (linked-list head, max coordinate).
// ChrNameMap: small flat-vector map keyed by chromosome name.
// Arena:      bump allocator for line tails / weight strings on stdin/gzip.
// Constants:  weight-field cap, bucket-sort default budget, hybrid cutoff.
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
	int len;			// max value of the bucket-sort key seen on this chromosome:
						//   beg for --sort s/b, or strand-aware 5'-end for --sort 5
	uint32_t lastRead;	// head of the linked list of read indices for this chromosome
						// (0 = no reads yet / end-of-list sentinel)
	uint32_t idx;		// sequential chromosome index, assigned after parsing & alphabetical sort
	chrInfoT() : len(0), lastRead(0), idx(0) {}
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

// Stack buffer used for the BED weight / RAL weight field (field 4 in BED,
// field 5 in RAL). Real weights are short numeric strings; 256 B is generous.
// On overflow, copyField returns 0 and the parser errors with a clear message
// rather than silently truncating.
const int kWeightBufSize = 256;

// Default per-allocation budget for the bucket-sort path's chromTable when
// `--max-mem` isn't set. 4 GB matches the previous fixed chrLenLimit of
// 1 Gbp × 4 B/slot — past this the bucket-sort allocation becomes a poor
// trade-off vs. the --low-mem-ssd path. Override via --max-mem=N[GMK].
static constexpr size_t kDefaultBucketChromBudget = 4ULL * 1024 * 1024 * 1024;

// Hybrid sort strategy cutoff (within the classic path): files with fewer
// than this many reads use index-array std::sort (O(n log n)). Larger files
// use bucket/counting sort (O(n + m) where m = max chromosome coordinate).
// Bucket sort allocates a chromTable[maxChrLen] array (up to 4 GB) on each
// thread, so it's RAM-hungry; the classic path is mostly used for small
// inputs anyway since --low-mem-ssd is faster and lighter at scale.
//
// The 200M default is conservative — std::sort beats bucket on the headline
// 200M-row fixture at -t 1 even though bucket has better asymptotic
// complexity. --bucket-cutoff at runtime lets users force either.
const int defaultBucketCutoff = 200000000;


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
// PARSERS — FIELDS, BED / RAL LINES, COLLAPSE-MODE WEIGHTS
// ----------------------------------------------------------------------------
// Field-level (skipField/parseUInt/copyField) and line-level
// (parseBedLine3/parseBedLineFull/parseRalLine) primitives, all `inline` for
// the hot path. Plus parseWeight/sumWeights* which the bucket-sort and
// classic-sort emit paths use to consume `--collapse` weight strings.
// Hand-written; faster than sscanf.
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

// RAL line parser. Records chr as pointer+length into the line buffer
// (caller owns the buffer; no copy / no fixed length limit) — same pattern
// parseBedLine3/parseBedLineFull use for chr.
static int parseRalLine(const char* buf,
                        const char** chrPtr, int* chrLen,
                        int* beg, int* end,
                        char* weightBuf, int weightMax,
                        char* strandChar)
{
	const char* p = buf;

	if (!skipField(&p)) return 0;

	// Field 2 (chr) — record pointer+length; no copying, no limit.
	*chrPtr = p;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	*chrLen = (int)(p - *chrPtr);
	if (*chrLen == 0) return 0;
	if (*p != '\t' && *p != ' ') return 0;
	p++;

	*strandChar = *p;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	if (*p != '\t' && *p != ' ') return 1;
	p++;

	int v = parseUInt(&p);
	if (v < 0) return 1;
	*beg = v;
	if (*p == '\t' || *p == ' ') p++;

	v = parseUInt(&p);
	if (v < 0) return 3;
	*end = v;

	if (*p != '\t' && *p != ' ') return 3;
	p++;
	int n = copyField(&p, weightBuf, weightMax);
	if (n == 0) return 3;
	if (n < 0) {
		fprintf(stderr, "Error: RAL weight field exceeds %d bytes in: %s\n",
		        weightMax - 1, buf);
		exit(1);
	}

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

// Sum weights across a contiguous buffer of read indices.
static inline float sumWeightsBuf(const seqread* reads, const uint32_t* idxs, int n)
{
	float sum = 0.0f;
	for(int j = 0; j < n; ++j) sum += parseWeight(reads, idxs[j]);
	return sum;
}

// Sum weights along a chromTable linked list at a given position, consuming the list.
// On return chromTable[pos] == 0.
static inline float sumWeightsList(seqread* reads, uint32_t* chromTable, int pos)
{
	float sum = 0.0f;
	while(chromTable[pos])
	{
		sum += parseWeight(reads, chromTable[pos]);
		chromTable[pos] = reads[chromTable[pos]].next;
	}
	return sum;
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
                          bool fRal, int fCollapse, char sortMode, int numThreads,
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
					char strandChar = '+';
					char weight[kWeightBufSize];
					weight[0] = '0'; weight[1] = '\0';
					const char* tailPtr = "";
					int numArgsRead;

					if(fRal)
					{
						numArgsRead = parseRalLine(linePtr, &chrPtr, &chrLen,
						                           &beg, &lineEnd, weight, kWeightBufSize,
						                           &strandChar);
					}
					else
					{
						numArgsRead = parseBedLine3(linePtr, &chrPtr, &chrLen,
						                            &beg, &lineEnd, &tailPtr);
					}
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
			char strandChar = '+';
			char weight[kWeightBufSize];
			weight[0] = '0'; weight[1] = '\0';
			const char* tailPtr = "";
			int numArgsRead;

			if(fRal)
			{
				numArgsRead = parseRalLine(linePtr, &chrPtr, &chrLen,
				                           &beg, &end, weight, kWeightBufSize,
				                           &strandChar);
			}
			else
			{
				numArgsRead = parseBedLine3(linePtr, &chrPtr, &chrLen, &beg, &end, &tailPtr);
			}

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
				int n;
				if(fRal)
					n = parseRalLine(linePtr, &chrPtr, &chrLen,
					                 &dummy_beg, &dummy_end, weight, kWeightBufSize, &dummy_str);
				else
					n = parseBedLineFull(linePtr, &chrPtr, &chrLen, &dummy_beg, &dummy_end,
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
				int n;
				if(fRal)
					n = parseRalLine(linePtr, &chrPtr, &chrLen,
					                 &dummy_beg, &dummy_end, weight, kWeightBufSize, &strandChar);
				else
					n = parseBedLineFull(linePtr, &chrPtr, &chrLen, &dummy_beg, &dummy_end,
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
// CLASSIC SORT PATH (default + bucket sort)
// ----------------------------------------------------------------------------
// Pipeline:
//   parse: parseLines<UseMmap> (serial / stdin)  OR  parseChunkMmap (parallel
//          via parseMmapDispatch / parseMmapSerial)
//   sort:  sortIndicesDispatch -> sortIndices<SortMode> -> radixSort64 or
//          std::sort with ReadCmp<SortMode>
//   emit:  per-chrom processChromBucketTpl<FullLine> for the bucket path,
//          or a flat scan over the sorted index for the classic path
//
// This section also owns the small write helpers (writeUInt, writeBedLine)
// used by both sort and bucket emit paths.
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
	// --sort b needs a tertiary key (end), which doesn't pack into 64 bits at
	// the BED coordinate range — falls through to std::sort. Same key-packing
	// works for both single-thread and multi-thread; radixSort64 internally
	// dispatches to a serial or parallel histogram + scatter pipeline.
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

// Bucket-sort one chromosome: scatter its linked list of reads into
// chromTable[chosenPosition], then scan positions in order and emit the result.
// chromTable and (for --sort b) numReadsBeg must be sized at least thisChrLen+1
// and zero-initialized; on return all touched slots are 0 again so the buffers
// can be reused for another chromosome.

// Templated on FullLine so gcc constant-folds the "full mmap line vs
// reconstruct from beg/end + tail" branch in the per-read emit loop;
// otherwise we'd pay an extra well-predicted-but-real fetch/decode branch
// per emitted line at -t 1 on multi-million-read files.
template<bool FullLine>
static __attribute__((always_inline)) inline void processChromBucketTpl(
    const char* chrName, int thisChrLen, uint32_t firstRead,
    seqread* reads, char sortMode, int fCollapse,
    uint32_t* chromTable, int* numReadsBeg, FILE* out)
{
	// Phase 1: scatter into chromTable linked lists keyed by chosenPosition.
	uint32_t currRead = firstRead;
	int maxNumReads = 0;
	while(currRead)
	{
		int chosenPosition = (sortMode == '5')
			? (reads[currRead].str == '-' ? reads[currRead].end : reads[currRead].beg)
			: reads[currRead].beg;
		uint32_t oldPtr = chromTable[chosenPosition];
		chromTable[chosenPosition] = currRead;
		currRead = reads[currRead].next;        // must read .next before overwriting it below
		reads[chromTable[chosenPosition]].next = oldPtr;
		if(sortMode == 'b')
		{
			numReadsBeg[chosenPosition]++;
			if(numReadsBeg[chosenPosition] > maxNumReads) maxNumReads++;
		}
	}

	// Phase 2: scan positions in order, emit and consume each bucket.
	if(sortMode == 'b')
	{
		uint32_t* readBuff = (uint32_t*) calloc((size_t)maxNumReads, sizeof(uint32_t));
		for(int pos = 0; pos <= thisChrLen; pos++)
		{
			if(!chromTable[pos]) continue;
			int n = 0;
			int foo = numReadsBeg[pos];
			while(chromTable[pos])
			{
				readBuff[n++] = chromTable[pos];
				chromTable[pos] = reads[chromTable[pos]].next;
			}
			std::sort(readBuff, readBuff + foo, [reads](uint32_t a, uint32_t b) {
				return reads[a].end < reads[b].end;
			});
			if(fCollapse)
				fprintf(out, "%s\t%d\t%d\t.\t%g\t+\n", chrName, pos, pos+1,
				        sumWeightsBuf(reads, readBuff, foo));
			else if constexpr(FullLine)
				for(int j = 0; j < foo; ++j)
				{
					fputs_unlocked(reads[readBuff[j]].line, out);
					fputc_unlocked('\n', out);
				}
			else
				for(int j = 0; j < foo; ++j)
					writeBedLine(chrName, reads[readBuff[j]].beg, reads[readBuff[j]].end,
					             reads[readBuff[j]].line, out);
			numReadsBeg[pos] = 0;
		}
		free(readBuff);
	}
	else  // sortMode in {5, s}: positions already give the final order.
	{
		for(int pos = 0; pos <= thisChrLen; pos++)
		{
			if(!chromTable[pos]) continue;
			if(fCollapse)
			{
				fprintf(out, "%s\t%d\t%d\t.\t%g\t+\n", chrName, pos, pos+1,
				        sumWeightsList(reads, chromTable, pos));
			}
			else if constexpr(FullLine)
			{
				while(chromTable[pos])
				{
					uint32_t ri = chromTable[pos];
					fputs_unlocked(reads[ri].line, out);
					fputc_unlocked('\n', out);
					chromTable[pos] = reads[ri].next;
				}
			}
			else
			{
				while(chromTable[pos])
				{
					uint32_t ri = chromTable[pos];
					writeBedLine(chrName, reads[ri].beg, reads[ri].end, reads[ri].line, out);
					chromTable[pos] = reads[ri].next;
				}
			}
		}
	}
}

// Runtime dispatch over the FullLine template parameter.
static __attribute__((always_inline)) inline void processChromBucket(
    const char* chrName, int thisChrLen, uint32_t firstRead,
    seqread* reads, char sortMode, int fCollapse, bool fullLine,
    uint32_t* chromTable, int* numReadsBeg, FILE* out)
{
	if(fullLine)
		processChromBucketTpl<true>(chrName, thisChrLen, firstRead, reads, sortMode, fCollapse,
		                            chromTable, numReadsBeg, out);
	else
		processChromBucketTpl<false>(chrName, thisChrLen, firstRead, reads, sortMode, fCollapse,
		                             chromTable, numReadsBeg, out);
}

// Parsing loop, templated on UseMmap to eliminate per-line branch.
// For BED, stores only the "tail" (fields 4+) in reads[].line (mmap: full line pointer;
// stdin/gzip: only tail copied to arena, saving ~50% arena memory).
// For RAL, stores the full line (output format differs from BED).
template<bool UseMmap>
static void parseLines(
    char* mmapBase, size_t mmapSize,
    FILE* fh,
    seqread*& reads, size_t& readCount, size_t& currMaxReads,
    string2chrInfoT& chrInfo, string2chrInfoT::iterator& thisChrIt,
    Arena* arena,
    bool fRal, int fCollapse, char sortMode)
{
	char* mmapCur = mmapBase;
	char* mmapLim = mmapBase + mmapSize;
	// getline() reuses and grows this buffer across calls, so there's no
	// fixed line-length limit on the stdin/gzip path. Only allocated when
	// UseMmap is false; freed at end of function.
	char* lineBuf = NULL;
	size_t lineBufCap = 0;
	const bool needExtra = fRal || fCollapse || (sortMode == '5');

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
		// reads[].next, the order[] array, chromTable[] entries, and radix-sort
		// internals are all uint32_t. Index 0 is reserved as the end-of-list
		// sentinel, so usable range is 1..UINT32_MAX-1.
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

		if(fRal)
		{
			numArgsRead = parseRalLine(linePtr, &chrPtr, &chrLen,
			                           &beg, &end, weight, kWeightBufSize,
			                           &strandChar);
		}
		else if(needExtra)
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
			if(fRal)
			{
				if(UseMmap)
					reads[readCount].line = linePtr;
				else
					reads[readCount].line = arena->alloc(linePtr, strlen(linePtr) + 1);
			}
			else if(!fCollapse)
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
			int chosenPosition = (sortMode == '5') ? (strandChar == '-' ? end : beg) : beg;

			// Fast same-chr check: compare length first, then bytes.
			if((int)thisChrIt->first.size() != chrLen ||
			   memcmp(thisChrIt->first.data(), chrPtr, chrLen) != 0)
			{
				// tryEmplaceByPtr avoids std::string construction on the common
				// "chrom already exists" path; only allocates on a true insert.
				auto ins = chrInfo.tryEmplaceByPtr(chrPtr, chrLen);
				thisChrIt = ins.first;
				if(!ins.second)
				{
					if(thisChrIt->second.len < chosenPosition) thisChrIt->second.len = chosenPosition;
					reads[readCount].next = thisChrIt->second.lastRead;
				}
				else
				{
					reads[readCount].next = 0;
					thisChrIt->second.len = chosenPosition;
				}
			}
			else
			{
				if(thisChrIt->second.len < chosenPosition) thisChrIt->second.len = chosenPosition;
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
	int len;         // max chosenPosition seen
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
                           bool fRal, int fCollapse, char sortMode,
                           Arena* arena)
{
	const bool needExtra = fRal || fCollapse || (sortMode == '5');
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

		if(fRal)
		{
			numArgsRead = parseRalLine(linePtr, &chrPtr, &chrLen,
			                           &beg, &lineEnd, weight, kWeightBufSize,
			                           &strandChar);
		}
		else if(needExtra)
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
		int chosenPosition = (sortMode == '5') ? (strandChar == '-' ? lineEnd : beg) : beg;

		// Same-chr fast path within this chunk.
		if(thisChrIt != result.chr.end() &&
		   (int)thisChrIt->first.size() == chrLen &&
		   memcmp(thisChrIt->first.data(), chrPtr, chrLen) == 0)
		{
			reads[idx].next = thisChrIt->second.head;
			thisChrIt->second.head = idx;
			if(thisChrIt->second.len < chosenPosition) thisChrIt->second.len = chosenPosition;
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
				thisChrIt->second.len = chosenPosition;
			}
			else
			{
				reads[idx].next = thisChrIt->second.head;
				thisChrIt->second.head = idx;
				if(thisChrIt->second.len < chosenPosition) thisChrIt->second.len = chosenPosition;
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
                                bool fRal, int fCollapse, char sortMode,
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
	                 fRal, fCollapse, sortMode);
	chrInfo.erase("");
	outReadCount = readCount - 1;
	return reads;
}

// Top-level mmap parser. Emits any leading header lines to stdout, then
// either runs parseMmapSerial (N==1) or chunks the body and parses in
// parallel. Returns the malloc'd reads array; caller frees.
static seqread* parseMmapDispatch(char* mmapBase, size_t mmapSize,
                                  bool fRal, int fCollapse, char sortMode,
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
		                       fRal, fCollapse, sortMode, arena,
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

	// reads[].next, the order[] array, chromTable[] entries, and radix-sort
	// internals are all uint32_t. Index 0 is reserved as the end-of-list
	// sentinel, so usable range is 1..UINT32_MAX-1.
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
			               fRal, fCollapse, sortMode, a);
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
				info.len = partial.len;
			}
			else
			{
				reads[partial.tail].next = info.lastRead;
				info.lastRead = partial.head;
				if(info.len < partial.len) info.len = partial.len;
			}
		}
	}

	outReadCount = totalReads;
	return reads;
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
		"Uses one thread by default and sorts at ~disk IO throughput limits.\n\n"
		"Input file should contain a new line character in the end of the last line.\n"
		"Gzip-compressed input (.gz) is transparently decompressed via gzip.\n"
		"BED header lines (track/browser/#) are passed through unchanged.\n\n"
		"Compilation-time limits:\n"
		"  Line length: no fixed limit (stdin/gzip uses getline; mmap uses memchr).\n"
		"  Chromosome name: no fixed limit (stored as pointer+length).\n"
		"  Weight field (BED col 4 / RAL): " + to_string(kWeightBufSize - 1) + " B; over-long values are rejected with an error.\n"
		"  Bucket-sort path: rejects a chromosome whose chromTable would exceed\n"
		"    --max-mem (default " + to_string(kDefaultBucketChromBudget >> 30) + " GB); use --low-mem-ssd otherwise.\n"
		"  Read number limit: " + to_string((unsigned long long)UINT32_MAX - 1) + " (all paths).\n"};
	app.set_version_flag("-V,--version", VERSION_STRING);

	string inputFile;
	char sortMode = 's';
	int fCollapse = 0;
	int fRal = 0;
	int lowMemSSD = 0;
	int bucketCutoff = defaultBucketCutoff;
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
	app.add_flag("-r,--ral", fRal, "input is in RAL format");
	app.add_flag("--collapse", fCollapse,
		"collapse BEDWEIGHT by summing weights of the reads and truncate the coordinates, "
		"set strand to \"+\", id to \".\". Makes no sense for --sort=b");
	app.add_flag("--low-mem-ssd", lowMemSSD,
		"low-memory two-pass file mode (SSD-friendly): keeps line offsets in RAM and "
		"sorts one chromosome at a time; slower than default but uses less peak RAM");
	app.add_option("--bucket-cutoff", bucketCutoff,
		"use bucket sort for files with at least this many reads; "
		"smaller files use std::sort (0 = always bucket sort)")
		->default_val(defaultBucketCutoff);
	app.add_option("-t,--threads", numThreads,
		"number of threads for the classic sort path "
		"(0 = all cores; 1 = single-threaded std::sort)")
		->default_val(0);
	app.add_flag("-n,--natural-sort", naturalSort,
		"sort chromosomes in natural order (chr2 < chr10) instead of "
		"lexicographic order (chr10 < chr2)");
	app.add_option("--max-mem", maxMemStr,
		"memory budget for the parallel bucket-sort path "
		"(e.g. 4G, 500M, 2048K, or bare bytes). Caps concurrent per-chromosome "
		"chromTable allocations so peak RAM stays within budget. "
		"Default: no cap (each thread allocates its own chromTable).");
	app.add_flag("-v,--verbose", verbose,
		"print parsing / sorting timing and per-chromosome length info to stderr "
		"(default: silent)");

	CLI11_PARSE(app, argc, argv);
	const size_t maxMemBytes = parseMemSize(maxMemStr);

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

	if(numThreads <= 0)
	{
		unsigned hw = std::thread::hardware_concurrency();
		numThreads = hw ? (int)hw : 1;
	}
	// TBB owns the worker pool used by std::execution::par. Capping it here
	// gives `--threads N` the same behavior the OpenMP build had.
	tbb::global_control tbbThreadCap(
		tbb::global_control::max_allowed_parallelism, (size_t)numThreads);

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
		if(!useMmap)
		{
			cerr << "Error: --low-mem-ssd requires file input (not stdin or gzip)" << endl;
			if(isGzip && fh) pclose(fh);
			return 1;
		}
		return lowMemSortMmap(mmapBase, mmapSize, fRal, fCollapse, sortMode,
		                      numThreads, naturalSort, maxMemBytes, verbose);
	}

	time_t tstart, tend;
	if(verbose) time(&tstart);

	int maxChrLen = 0;
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
		                          (bool)fRal, fCollapse, sortMode,
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
		                  fRal, fCollapse, sortMode);
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
		if(verbose) cerr << it->first << ": " << it->second.len << endl;
		chroms.push_back(it->first);
		if(maxChrLen < it->second.len)
			maxChrLen = it->second.len;
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
	bool useBucketSort = (bucketCutoff == 0) || (totalReads >= (size_t)bucketCutoff);
	if(verbose)
	{
		cerr << "We have " << totalReads << " regions, "
		     << chroms.size() << " chromosomes. "
		     << (useBucketSort ? "Using bucket sort." : "Using classic sort.")
		     << endl;
		time(&tstart);
	}

	if(!useBucketSort)
	{
		/*************************************************************
		 * CLASSIC SORT PATH — O(n log n) comparison-based sort
		 *
		 * Sort an index array (4-byte ints) rather than 24-byte seqread structs.
		 * Steps:
		 *  1. Walk linked lists, stamp chrIdx (reusing next's storage), build index array.
		 *  2. std::sort index array with inlined lambda comparator.
		 *  3. Linear scan to print (and optionally collapse).
		 *************************************************************/

		uint32_t* order = (uint32_t*) malloc(totalReads * sizeof(uint32_t));
		size_t oi = 0;
		for(uint32_t ci = 0; ci < (uint32_t)chroms.size(); ci++)
		{
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

		// Both 1-thread and N-thread paths go through sortIndicesDispatch which
		// picks the templated comparator at compile time and chooses serial vs
		// std::execution::par based on numThreads.
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
		else
		{
			// Hoist fullLine outside the per-read loop so gcc fully constant-folds it.
			if(fRal || useMmap)
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
		}
		free(order);
	}
	else
	{
		/*************************************************************
		 * BUCKET SORT PATH — O(n + m) counting/bucket sort
		 *
		 * Single-thread: one shared chromTable[maxChrLen+1] reused across chromosomes.
		 * Multi-thread: chromosomes processed in parallel via std::for_each(par, ...).
		 *   Each worker allocates its own per-chrom chromTable, writes output to an
		 *   open_memstream buffer, then waits at a producer-consumer barrier for its
		 *   alphabetical turn before flushing to stdout. The chromTable is freed
		 *   BEFORE the wait so a thread holding chr1 doesn't pin ~1 GB while it waits.
		 *************************************************************/
		const bool fullLine = fRal || useMmap;

		if(numThreads <= 1)
		{
			// One chromTable shared across all chroms, sized to the largest.
			// Reject before allocating if the slab would exceed the budget
			// (--max-mem if set, else 4 GB). --low-mem-ssd handles arbitrary
			// chromosome lengths since it never builds a per-position slab.
			size_t slabBytes = ((size_t)maxChrLen + 1) * sizeof(int);
			if(sortMode == 'b') slabBytes *= 2;
			size_t bucketBudget = (maxMemBytes > 0) ? maxMemBytes : kDefaultBucketChromBudget;
			if(slabBytes > bucketBudget)
			{
				cerr << "Error: largest chromosome has positions up to " << maxChrLen
				     << ", which would require " << (slabBytes >> 30) << " GB for the "
				     << "bucket-sort slab (budget: " << (bucketBudget >> 30) << " GB). "
				     << "Re-run with --low-mem-ssd, or raise --max-mem." << endl;
				exit(1);
			}
			uint32_t* chromTable = (uint32_t*) calloc((size_t)maxChrLen + 1, sizeof(uint32_t));
			int* numReadsBeg = (sortMode == 'b')
				? (int*) calloc((size_t)maxChrLen + 1, sizeof(int))
				: NULL;
			for(const std::string& chrom : chroms)
			{
				string2chrInfoT::iterator cit = chrInfo.find(chrom);
				processChromBucket(chrom.c_str(), cit->second.len, cit->second.lastRead,
				                   reads, sortMode, fCollapse, fullLine,
				                   chromTable, numReadsBeg, stdout);
			}
			free(chromTable);
			if(numReadsBeg) free(numReadsBeg);
		}
		else
		{
			std::vector<int> idxs(chroms.size());
			for(size_t i = 0; i < chroms.size(); i++) idxs[i] = (int)i;

			std::mutex printMtx;
			std::condition_variable printCv;
			int nextChromToPrint = 0;

			// Optional chromTable memory-budget gate. When --max-mem is set,
			// each task acquires `effCost` bytes from `budgetRemaining` before
			// allocating its chromTable, and releases them after the chrom is
			// fully processed. If a single chrom's cost exceeds the budget it
			// runs alone (clamping its draw to the full budget).
			std::mutex memMtx;
			std::condition_variable memCv;
			size_t budgetRemaining = maxMemBytes;

			std::for_each(std::execution::par, idxs.begin(), idxs.end(),
				[&](int ci) {
					const std::string& chrom = chroms[ci];
					string2chrInfoT::iterator cit = chrInfo.find(chrom);
					int thisChrLen = cit->second.len;

					// Cost = chromTable + (numReadsBeg if --sort b).
					size_t cost = (size_t)(thisChrLen + 1) * sizeof(int);
					if(sortMode == 'b') cost *= 2;

					// Reject before allocating if a single chrom's slab exceeds
					// the budget. --max-mem if set is the budget; otherwise
					// fall back to the 4 GB default. --low-mem-ssd handles
					// arbitrary chromosome lengths.
					size_t bucketBudget = (maxMemBytes > 0) ? maxMemBytes : kDefaultBucketChromBudget;
					if(cost > bucketBudget)
					{
						cerr << "Error: chromosome " << chrom << " has positions up to "
						     << thisChrLen << ", which would require " << (cost >> 30)
						     << " GB for the bucket-sort slab (budget: "
						     << (bucketBudget >> 30) << " GB). Re-run with --low-mem-ssd, "
						     << "or raise --max-mem." << endl;
						exit(1);
					}
					size_t effCost = (maxMemBytes > 0)
						? std::min(cost, maxMemBytes)  // a chrom bigger than the
						                               // whole budget runs alone
						: 0;

					if(effCost > 0)
					{
						std::unique_lock<std::mutex> lk(memMtx);
						memCv.wait(lk, [&]{ return budgetRemaining >= effCost; });
						budgetRemaining -= effCost;
					}

					char* memBuf = NULL;
					size_t memBufSz = 0;
					FILE* sink = open_memstream(&memBuf, &memBufSz);
					if(!sink)
					{
						cerr << "open_memstream failed for chrom " << chrom << endl;
						exit(1);
					}
					{
						std::vector<uint32_t> chromTable((size_t)thisChrLen + 1, 0);
						std::vector<int> numReadsBeg;
						if(sortMode == 'b') numReadsBeg.assign((size_t)thisChrLen + 1, 0);

						processChromBucket(chrom.c_str(), thisChrLen, cit->second.lastRead,
						                   reads, sortMode, fCollapse, fullLine,
						                   chromTable.data(),
						                   sortMode == 'b' ? numReadsBeg.data() : NULL,
						                   sink);
						// chromTable / numReadsBeg destruct here, releasing ~thisChrLen*4
						// bytes BEFORE we block at the print barrier.
					}
					fclose(sink);

					if(effCost > 0)
					{
						{
							std::lock_guard<std::mutex> lk(memMtx);
							budgetRemaining += effCost;
						}
						memCv.notify_all();
					}

					{
						std::unique_lock<std::mutex> lk(printMtx);
						printCv.wait(lk, [&]{ return nextChromToPrint == ci; });
						fwrite_unlocked(memBuf, 1, memBufSz, stdout);
						free(memBuf);
						nextChromToPrint++;
					}
					printCv.notify_all();
				});
		}
	}
	if(verbose)
	{
		time(&tend);
		fprintf(stderr, "Sorting has taken %d seconds\n",
		        (int)(difftime(tend, tstart)+0.5));
	}
	return 0;
}
