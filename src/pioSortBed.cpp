#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <algorithm>
#include <condition_variable>
#include <execution>
#include <iostream>
#include <mutex>
#include <thread>
#include <time.h>
#include <unordered_map>
#include <tbb/global_control.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "CLI11.hpp"

#ifndef VERSION_STRING
#define VERSION_STRING "2.0.0"
#endif

/* Copyright: Piotr Balwierz */

using namespace std;

// Field layout minimizes padding: pointer first, then ints, then char.
// 24 bytes per read on 64-bit (8 + 4 + 4 + 4 + 1 + 3 padding).
//
// The `next` field (linked lists during parsing) and `chrIdx` field (classic
// sort path integer chromosome index) share the same 4 bytes via a union —
// the two sort paths never need both simultaneously.
class seqread
{
	public:
	char* line; // keeps the full line; in case of collapse option just the weight as string.
	int beg;
	int end;
	union { int next; int chrIdx; };
	char str;
};

// make sure we call find or [] on a chr only once!
class chrInfoT
{
	public:
	int len;		// length (the max beg position)
	int lastRead;	// keeps the head of the list of all reads at a given chromosome
	int idx;		// sequential chromosome index, assigned after parsing & alphabetical sort
	chrInfoT() : len(0), lastRead(0), idx(0) {}
};

typedef std::unordered_map<string, chrInfoT> string2chrInfoT;

const int lineBufSize = 1024;		// buffer for the read file lines (stdin path only; mmap has no line length limit)
const int chrLenLimit = 1000000000;	// 1 Gbp. we don't believe there are longer chromosomes;
									// Just a safety check, not to allocate TBytes of ram
									// Increase if you know what you are doing.
									// And change the text chrTooLongMsg below.
									// Corresponds to 1000M * 4B = 4GB on a 64bit machine
const char chrTooLongMsg[] = "That is more than 3 times the length of human chr1! If you are sure this is correct recompile the program with an increased chrLenLimit value";
const int chrNameBufSize = 256;		// buffer for chromosme names like "chr1_random"
									// and strand too, although it should be only "+" or "-"

// Hybrid sort strategy cutoff: files with fewer than this many reads use
// classic O(n log n) comparison sort (std::sort on an index array).
// Larger files use bucket/counting sort which is O(n + m) where m = max
// chromosome length.  Bucket sort allocates a chromTable[maxChrLen] array
// (up to 4 GB), so the classic path saves memory for smaller files while
// bucket sort wins on throughput for large genomic datasets.
const int defaultBucketCutoff = 50000000;


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

static inline int copyField(const char** pp, char* dst, int dstMax)
{
	const char* p = *pp;
	int n = 0;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n')
	{
		if (n < dstMax - 1) dst[n++] = *p;
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

// Full BED parser: chr, beg, end, weight (field 4), strand (field 5).
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

	// field 3: id — skip it
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	if (*p != '\t' && *p != ' ') return 3;
	p++;

	// field 4: weight
	int n = copyField(&p, weightBuf, weightMax);
	if (n == 0) return 3;

	// field 5: strand — copyField already consumed the tab before this field
	if (!*p || *p == '\r' || *p == '\n') return 4;
	*strandChar = *p;
	return 5;
}

// RAL line parser.
static int parseRalLine(const char* buf,
                        char* chr, int chrMax,
                        int* beg, int* end,
                        char* weightBuf, int weightMax,
                        char* strandChar)
{
	const char* p = buf;

	if (!skipField(&p)) return 0;

	int n = copyField(&p, chr, chrMax);
	if (n == 0) return 0;

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
	n = copyField(&p, weightBuf, weightMax);
	if (n == 0) return 3;

	return 5;
}

// Parse a weight string written into reads[idx].line (collapse mode) using
// strtof (faster + locale-independent compared to sscanf). Aborts on malformed input.
static inline float parseWeight(const seqread* reads, int idx)
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

// Sum weights across a contiguous int buffer of read indices.
static inline float sumWeightsBuf(const seqread* reads, const int* idxs, int n)
{
	float sum = 0.0f;
	for(int j = 0; j < n; ++j) sum += parseWeight(reads, idxs[j]);
	return sum;
}

// Sum weights along a chromTable linked list at a given position, consuming the list.
// On return chromTable[pos] == 0.
static inline float sumWeightsList(seqread* reads, int* chromTable, int pos)
{
	float sum = 0.0f;
	while(chromTable[pos])
	{
		sum += parseWeight(reads, chromTable[pos]);
		chromTable[pos] = reads[chromTable[pos]].next;
	}
	return sum;
}

struct lowMemNode
{
	size_t off;
	int next;
};

struct lowMemRec
{
	int beg;
	int end;
	char str;
	float weight;
	const char* line;
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

// Low-memory file mode (SSD-friendly):
// Pass 1: walk mmap once, build per-chromosome linked lists of line offsets.
// Pass 2: process one chromosome at a time (parse, sort, print), so peak RAM
// depends on the largest chromosome chunk, not the whole file.
static int lowMemSortMmap(char* mmapBase, size_t mmapSize,
                          bool fRal, int fCollapse, char sortMode, int numThreads,
                          bool naturalSort)
{
	time_t tstart, tend;
	time(&tstart);

	string2chrInfoT chrInfo;
	pair<string2chrInfoT::iterator,bool> insResult =
		chrInfo.insert(make_pair(string(""), chrInfoT()));
	string2chrInfoT::iterator thisChrIt = insResult.first;

	int nodeCap = 1024;
	int nodeCount = 1; // index 0 reserved as end-of-list sentinel
	lowMemNode* nodes = (lowMemNode*) malloc(nodeCap * sizeof(lowMemNode));
	if(!nodes)
	{
		cerr << "Error: out of memory allocating low-memory offset table" << endl;
		return 1;
	}

	// Pass 1: index all lines by chromosome.
	char* mmapCur = mmapBase;
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

		// Pass through BED header lines (track/browser/# comments) immediately.
		if(linePtr[0] == '#' ||
		   strncmp(linePtr, "track ", 6) == 0 ||
		   strncmp(linePtr, "browser ", 8) == 0)
		{
			fputs_unlocked(linePtr, stdout);
			fputc_unlocked('\n', stdout);
			continue;
		}

		const char* chrPtr;
		int chrLen;
		int beg = 0;
		int end = 0;
		int numArgsRead;
		char strandChar = '+';
		char weight[chrNameBufSize];
		weight[0] = '0'; weight[1] = '\0';
		const char* tailPtr = "";
		char chrBuf[chrNameBufSize];

		if(fRal)
		{
			numArgsRead = parseRalLine(linePtr, chrBuf, chrNameBufSize,
			                           &beg, &end, weight, chrNameBufSize,
			                           &strandChar);
			chrPtr = chrBuf;
			chrLen = (int)strlen(chrBuf);
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
			pair<string2chrInfoT::iterator,bool> ins =
				chrInfo.insert(make_pair(string(chrPtr, chrLen), chrInfoT()));
			thisChrIt = ins.first;
		}

		if(nodeCount == nodeCap)
		{
			nodeCap = (int)(nodeCap * 2);
			lowMemNode* tmp = (lowMemNode*) realloc(nodes, nodeCap * sizeof(lowMemNode));
			if(!tmp)
			{
				cerr << "Error: out of memory growing low-memory offset table" << endl;
				free(nodes);
				return 1;
			}
			nodes = tmp;
		}

		nodes[nodeCount].off = (size_t)(linePtr - mmapBase);
		nodes[nodeCount].next = thisChrIt->second.lastRead;
		thisChrIt->second.lastRead = nodeCount;
		nodeCount++;
	}

	chrInfo.erase("");
	vector<string> chroms;
	chroms.reserve(chrInfo.size());
	for(string2chrInfoT::iterator it = chrInfo.begin(); it != chrInfo.end(); it++)
		chroms.push_back(it->first);
	if(naturalSort)
		std::sort(chroms.begin(), chroms.end(), naturalChrLess);
	else
		std::sort(chroms.begin(), chroms.end());

	time(&tend);
	fprintf(stderr, "Reading has taken %d seconds\n", (int)(difftime(tend, tstart)+0.5));
	cerr << "We have " << nodeCount - 1 << " regions.\n"
		 << chroms.size() << " chromosomes\nSorting..." << endl;
	time(&tstart);

	// Pass 2 accesses lines in chromosome-sorted order (not file order).
	madvise(mmapBase, mmapSize, MADV_RANDOM);

	// Pass 2: for each chromosome, re-parse lines, sort, and print.
	for(size_t ci = 0; ci < chroms.size(); ci++)
	{
		string2chrInfoT::iterator cit = chrInfo.find(chroms[ci]);
		vector<int> nodeIdx;
		for(int cur = cit->second.lastRead; cur; cur = nodes[cur].next)
			nodeIdx.push_back(cur);

		vector<lowMemRec> recs;
		recs.reserve(nodeIdx.size());
		for(size_t i = 0; i < nodeIdx.size(); i++)
		{
			const char* linePtr = mmapBase + nodes[nodeIdx[i]].off;
			int beg = 0, end = 0;
			char strandChar = '+';
			char weight[chrNameBufSize];
			weight[0] = '0'; weight[1] = '\0';
			const char* chrPtr;
			int chrLen;
			int numArgsRead;
			const char* tailPtr = "";
			char chrBuf[chrNameBufSize];

			if(fRal)
			{
				numArgsRead = parseRalLine(linePtr, chrBuf, chrNameBufSize,
				                           &beg, &end, weight, chrNameBufSize,
				                           &strandChar);
			}
			else if(fCollapse || sortMode == '5')
			{
				numArgsRead = parseBedLineFull(linePtr, &chrPtr, &chrLen,
				                              &beg, &end, weight, chrNameBufSize,
				                              &strandChar, &tailPtr);
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

			lowMemRec r;
			r.beg = beg;
			r.end = end;
			r.str = strandChar;
			r.line = linePtr;
			r.weight = 0.0f;
			if(fCollapse)
			{
				char* endptr;
				r.weight = strtof(weight, &endptr);
				if(endptr == weight)
				{
					cerr << "Malformed weight " << weight << endl;
					free(nodes);
					return 1;
				}
			}
			recs.push_back(r);
		}

		auto cmp = [sortMode](const lowMemRec& a, const lowMemRec& b) {
			if(sortMode == 'b')
			{
				if(a.beg != b.beg) return a.beg < b.beg;
				return a.end < b.end;
			}
			if(sortMode == '5')
			{
				int pa = (a.str == '-') ? a.end : a.beg;
				int pb = (b.str == '-') ? b.end : b.beg;
				return pa < pb;
			}
			return a.beg < b.beg;
		};

		if(numThreads == 1)
			std::sort(recs.begin(), recs.end(), cmp);
		else
			std::sort(std::execution::par, recs.begin(), recs.end(), cmp);

		if(fCollapse)
		{
			size_t i = 0;
			while(i < recs.size())
			{
				int pos = recs[i].beg;
				float sum = 0.0f;
				while(i < recs.size() && recs[i].beg == pos)
				{
					sum += recs[i].weight;
					i++;
				}
				printf("%s\t%d\t%d\t.\t%g\t+\n", chroms[ci].c_str(), pos, pos+1, sum);
			}
		}
		else
		{
			for(size_t i = 0; i < recs.size(); i++)
			{
				fputs_unlocked(recs[i].line, stdout);
				fputc_unlocked('\n', stdout);
			}
		}
	}

	free(nodes);
	time(&tend);
	fprintf(stderr, "Sorting has taken %d seconds\n", (int)(difftime(tend, tstart)+0.5));
	return 0;
}

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

// Print one read. fullLine==true means reads[idx].line holds the entire input
// line (RAL or mmap path) — a single fputs+\n is enough. fullLine==false means
// .line is just the tail (stdin/gzip BED path), so reconstruct chr/beg/end.
static inline void printRead(const seqread* reads, int idx, const char* chr,
                             bool fullLine, FILE* out)
{
	if(fullLine)
	{
		fputs_unlocked(reads[idx].line, out);
		fputc_unlocked('\n', out);
	}
	else
	{
		writeBedLine(chr, reads[idx].beg, reads[idx].end, reads[idx].line, out);
	}
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
	bool operator()(int a, int b) const
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

template<char SortMode>
static void sortIndices(int* order, int n, seqread* reads, int numThreads)
{
	ReadCmp<SortMode> cmp{reads};
	if(numThreads == 1)
		std::sort(order, order + n, cmp);
	else
		std::sort(std::execution::par, order, order + n, cmp);
}

// Dispatch the templated sort by the runtime sortMode char.
static void sortIndicesDispatch(int* order, int n, seqread* reads, char sortMode, int numThreads)
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
static void processChromBucket(const char* chrName, int thisChrLen, int firstRead,
                               seqread* reads, char sortMode, int fCollapse, bool fullLine,
                               int* chromTable, int* numReadsBeg, FILE* out)
{
	// Phase 1: scatter into chromTable linked lists keyed by chosenPosition.
	int currRead = firstRead;
	int maxNumReads = 0;
	while(currRead)
	{
		int chosenPosition = (sortMode == '5')
			? (reads[currRead].str == '-' ? reads[currRead].end : reads[currRead].beg)
			: reads[currRead].beg;
		int oldPtr = chromTable[chosenPosition];
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
		int* readBuff = (int*) calloc((size_t)maxNumReads, sizeof(int));
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
			std::sort(readBuff, readBuff + foo, [reads](int a, int b) {
				return reads[a].end < reads[b].end;
			});
			if(fCollapse)
				fprintf(out, "%s\t%d\t%d\t.\t%g\t+\n", chrName, pos, pos+1,
				        sumWeightsBuf(reads, readBuff, foo));
			else
				for(int j = 0; j < foo; ++j)
					printRead(reads, readBuff[j], chrName, fullLine, out);
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
				fprintf(out, "%s\t%d\t%d\t.\t%g\t+\n", chrName, pos, pos+1,
				        sumWeightsList(reads, chromTable, pos));
			else
				while(chromTable[pos])
				{
					int ri = chromTable[pos];
					printRead(reads, ri, chrName, fullLine, out);
					chromTable[pos] = reads[ri].next;
				}
		}
	}
}

// Parsing loop, templated on UseMmap to eliminate per-line branch.
// For BED, stores only the "tail" (fields 4+) in reads[].line (mmap: full line pointer;
// stdin/gzip: only tail copied to arena, saving ~50% arena memory).
// For RAL, stores the full line (output format differs from BED).
template<bool UseMmap>
static void parseLines(
    char* mmapBase, size_t mmapSize,
    FILE* fh,
    seqread*& reads, int& readCount, int& currMaxReads,
    string2chrInfoT& chrInfo, string2chrInfoT::iterator& thisChrIt,
    Arena* arena,
    bool fRal, int fCollapse, char sortMode)
{
	char* mmapCur = mmapBase;
	char* mmapLim = mmapBase + mmapSize;
	char stdinBuf[lineBufSize];
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
			if(!fgets_unlocked(stdinBuf, lineBufSize, fh)) break;
			linePtr = stdinBuf;
			size_t len = strlen(linePtr);
			if (len >= 2 && linePtr[len-2] == '\r' && linePtr[len-1] == '\n')
				{ linePtr[len-2] = '\0'; }
			else if (len >= 1 && linePtr[len-1] == '\n')
				{ linePtr[len-1] = '\0'; }
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
			currMaxReads = (int) (currMaxReads * 2);
			seqread* tmp = (seqread*) realloc(reads, currMaxReads * sizeof(seqread));
			if(!tmp) { fprintf(stderr, "Error: out of memory\n"); exit(1); }
			reads = tmp;
		}

		int beg = 0;
		int end = 0;
		char strandChar = '+';
		char weight[chrNameBufSize];
		weight[0] = '0'; weight[1] = '\0';
		char chrBuf[chrNameBufSize]; // only used by RAL parser
		const char* chrPtr;
		int chrLen;
		const char* tailPtr = "";
		int numArgsRead;

		if(fRal)
		{
			numArgsRead = parseRalLine(linePtr, chrBuf, chrNameBufSize,
			                           &beg, &end, weight, chrNameBufSize,
			                           &strandChar);
			chrPtr = chrBuf;
			chrLen = (int)strlen(chrBuf);
		}
		else if(needExtra)
		{
			numArgsRead = parseBedLineFull(linePtr, &chrPtr, &chrLen,
			                              &beg, &end, weight, chrNameBufSize,
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
				pair<string2chrInfoT::iterator,bool> ins =
					chrInfo.insert(make_pair(string(chrPtr, chrLen), chrInfoT()));
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
			thisChrIt->second.lastRead = readCount;
		}
		else
		{
			cerr << "Error in parsing line: " << linePtr << endl
				<< "Perhaps this line is malformed?" << endl;
			exit(1);
		}
		readCount ++;
	}
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
	int head;   // global idx of head (most recently inserted read)
	int tail;   // global idx of tail (read whose .next is 0)
	int len;    // max chosenPosition seen
};

typedef std::unordered_map<std::string, ChunkChrPartial> ChunkChrMap;

struct ChunkResult
{
	int firstIdx;
	int count;
	ChunkChrMap chr;
};

// Parse [start, end) of the mmap into reads[base+1 .. base+count].
// All indices stored (head/tail and reads[i].next) are global, so no
// rebase pass is needed. Mutates the chunk in place: replaces '\n' (and
// preceding '\r') with '\0' to NUL-terminate each line.
static void parseChunkMmap(char* start, char* end,
                           seqread* reads, int base,
                           ChunkResult& result,
                           bool fRal, char sortMode)
{
	const bool needExtra = fRal || (sortMode == '5');
	char* mmapCur = start;
	char* mmapLim = end;
	int idx = base + 1;
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
		char weight[chrNameBufSize];
		weight[0] = '0'; weight[1] = '\0';
		char chrBuf[chrNameBufSize];
		const char* chrPtr;
		int chrLen;
		const char* tailPtr = "";
		int numArgsRead;

		if(fRal)
		{
			numArgsRead = parseRalLine(linePtr, chrBuf, chrNameBufSize,
			                           &beg, &lineEnd, weight, chrNameBufSize,
			                           &strandChar);
			chrPtr = chrBuf;
			chrLen = (int)strlen(chrBuf);
		}
		else if(needExtra)
		{
			numArgsRead = parseBedLineFull(linePtr, &chrPtr, &chrLen,
			                               &beg, &lineEnd, weight, chrNameBufSize,
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

		reads[idx].line = linePtr;
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
			auto ins = result.chr.try_emplace(std::string(chrPtr, chrLen));
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
                                int& outReadCount)
{
	int currMaxReads = 1024;
	size_t bodySize = mmapSize - bodyStart;
	size_t scanSize = bodySize < 65536 ? bodySize : 65536;
	int nLines = 0;
	for(size_t i = 0; i < scanSize; i++)
		if(mmapBase[bodyStart + i] == '\n') nLines++;
	if(nLines > 0)
		currMaxReads = (int)((double)bodySize / (double)scanSize * nLines * 1.1);
	if(currMaxReads < 1024) currMaxReads = 1024;

	seqread* reads = (seqread*) malloc((size_t)currMaxReads * sizeof(seqread));
	if(!reads)
	{
		cerr << "Error: out of memory allocating reads array" << endl;
		exit(1);
	}

	int readCount = 1;  // slot 0 is the end-of-list sentinel
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
                                  string2chrInfoT& chrInfo,
                                  int& outReadCount)
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

	// --collapse mode in parallel would need per-thread arenas for the weight
	// strings — fall back to the serial parser for now.
	if(fCollapse)
		return parseMmapSerial(mmapBase, bodyStart, mmapSize,
		                       fRal, fCollapse, sortMode, arena,
		                       chrInfo, outReadCount);

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
	std::vector<int> chunkCount(N, 0);
	std::vector<int> indices(N);
	for(int i = 0; i < N; i++) indices[i] = i;
	std::for_each(std::execution::par, indices.begin(), indices.end(),
		[&](int i) {
			size_t a = chunkStart[i], b = chunkStart[i + 1];
			int c = 0;
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

	std::vector<int> chunkBase(N + 1);
	chunkBase[0] = 0;
	for(int i = 0; i < N; i++) chunkBase[i + 1] = chunkBase[i] + chunkCount[i];
	int totalReads = chunkBase[N];

	seqread* reads = (seqread*) malloc(((size_t)totalReads + 1) * sizeof(seqread));
	if(!reads)
	{
		cerr << "Error: out of memory allocating reads array (" << totalReads << " entries)" << endl;
		exit(1);
	}

	std::vector<ChunkResult> chunkRes(N);
	std::for_each(std::execution::par, indices.begin(), indices.end(),
		[&](int i) {
			parseChunkMmap(mmapBase + chunkStart[i], mmapBase + chunkStart[i + 1],
			               reads, chunkBase[i], chunkRes[i],
			               fRal, sortMode);
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
		"  Line length limit: " + to_string(lineBufSize) + " (stdin only; no limit for file input)\n"
		"  Chromosome name limit: " + to_string(chrNameBufSize) + "\n"
		"  Chromosome length limit: " + to_string(chrLenLimit/1000000) + "Mbp\n"
		"  Chromosome/contig number limit: unlimited\n"
		"  Read number limit: unlimited (but will be kept in the memory)\n"};
	app.set_version_flag("-V,--version", VERSION_STRING);

	string inputFile;
	char sortMode = 's';
	int fCollapse = 0;
	int fRal = 0;
	int lowMemSSD = 0;
	int bucketCutoff = defaultBucketCutoff;
	int numThreads = 0;
	bool naturalSort = false;

	app.add_option("input-file", inputFile,
		"input file; \"-\" reads from stdin; .gz files are decompressed automatically")
		->required();
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

	CLI11_PARSE(app, argc, argv);

	if(numThreads <= 0)
	{
		unsigned hw = std::thread::hardware_concurrency();
		numThreads = hw ? (int)hw : 1;
	}
	// TBB owns the worker pool used by std::execution::par. Capping it here
	// gives `--threads N` the same behavior the OpenMP build had.
	tbb::global_control tbbThreadCap(
		tbb::global_control::max_allowed_parallelism, (size_t)numThreads);

	int readCount = 1;	// count == 0 is used as an end of a list.
	int currMaxReads = 1024;
	long fileSize = -1;

	char* mmapBase = NULL;
	size_t mmapSize = 0;
	bool useMmap = false;
	bool isGzip = false;
	FILE* fh = NULL;
	Arena* arena = NULL;

	if(inputFile == "-")
	{
		fh = stdin;
		cerr << "Reading data from standard input" << endl;
		const size_t chunkSize = 64UL * 1024 * 1024;
		arena = new Arena(chunkSize, chunkSize);
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
			cerr << "Reading gzip data from " << inputFile << endl;
			const size_t chunkSize = 64UL * 1024 * 1024;
			arena = new Arena(chunkSize, chunkSize);
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
			mmapBase = (char*) mmap(NULL, (size_t)fileSize + 1, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
			if(mmapBase == MAP_FAILED)
			{
				cerr << "Error: mmap failed for " << inputFile << endl;
				close(fd);
				return 1;
			}
			close(fd);
			madvise(mmapBase, (size_t)fileSize, MADV_SEQUENTIAL);
			mmapSize = (size_t)fileSize;
			useMmap = true;

			// Estimate read count from first 64KB
			size_t scanSize = mmapSize < 65536 ? mmapSize : 65536;
			int nLines = 0;
			for(size_t i = 0; i < scanSize; i++)
				if(mmapBase[i] == '\n') nLines++;
			if(nLines > 0)
				currMaxReads = (int)((double)mmapSize / scanSize * nLines * 1.1);
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
		                      numThreads, naturalSort);
	}

	time_t tstart, tend;
	time(&tstart);

	int maxChrLen = 0;
	string2chrInfoT chrInfo;
	seqread* reads = NULL;

	if(useMmap)
	{
		// Parallel parser dispatches to parseMmapSerial when N==1 or --collapse,
		// or to parseChunkMmap-per-chunk otherwise. Header lines (track/browser/#)
		// are emitted to stdout in a leading pre-pass before chunking.
		int parsed = 0;
		reads = parseMmapDispatch(mmapBase, mmapSize,
		                          (bool)fRal, fCollapse, sortMode,
		                          numThreads, arena, chrInfo, parsed);
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

	time(&tend);
	double difftime_reading = difftime(tend, tstart);
	fprintf(stderr, "Reading has taken %d seconds\n", (int)(difftime_reading+0.5));

	/***********************************
	 *  the actual sorting starts here *
	 * *********************************/
	int nChrom = chrInfo.size();
	std::vector<std::string> chroms;
	for(string2chrInfoT::iterator it=chrInfo.begin(); it!=chrInfo.end(); it++)
	{
		cerr << it->first << ": " << it->second.len << endl;
		chroms.push_back(it->first);
		if(maxChrLen < it->second.len)
			maxChrLen = it->second.len;
	}
	cerr << "We have " << readCount-1 << " regions.\n"
		<< nChrom << " chromosomes\nSorting..." << endl;
	if(maxChrLen > chrLenLimit)
	{
		cerr << "Error: There is a region starting at " << maxChrLen << "." << endl
			<<  chrTooLongMsg << endl;
		exit(1);
	}
	time(&tstart);

	if(naturalSort)
		std::sort(chroms.begin(), chroms.end(), naturalChrLess);
	else
		std::sort(chroms.begin(), chroms.end());

	// Assign sequential chromosome indices so the classic sort path can compare
	// chromosomes with a single int compare instead of a string compare.
	for(int ci = 0; ci < (int)chroms.size(); ci++)
		chrInfo.find(chroms[ci])->second.idx = ci;

	int totalReads = readCount - 1;
	bool useBucketSort = (bucketCutoff == 0) || (totalReads >= bucketCutoff);
	if(useBucketSort)
		cerr << "Using bucket sort (" << totalReads << " reads >= cutoff " << bucketCutoff << ")" << endl;
	else
		cerr << "Using classic sort (" << totalReads << " reads < cutoff " << bucketCutoff << ")" << endl;

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

		int* order = (int*) malloc(totalReads * sizeof(int));
		int oi = 0;
		for(int ci = 0; ci < (int)chroms.size(); ci++)
		{
			string2chrInfoT::iterator cit = chrInfo.find(chroms[ci]);
			int curr = cit->second.lastRead;
			while(curr)
			{
				int nxt = reads[curr].next;
				reads[curr].chrIdx = ci;
				order[oi++] = curr;
				curr = nxt;
			}
		}

		// Both 1-thread and N-thread paths go through sortIndicesDispatch which
		// picks the templated comparator at compile time and chooses serial vs
		// std::execution::par based on numThreads.
		if(numThreads > 1)
			cerr << "Parallel sort using " << numThreads << " threads" << endl;
		sortIndicesDispatch(order, totalReads, reads, sortMode, numThreads);

		if(fCollapse)
		{
			int i = 0;
			while(i < totalReads)
			{
				int ri = order[i];
				int ci = reads[ri].chrIdx;
				int pos = reads[ri].beg;
				float sum = 0.0f;
				while(i < totalReads)
				{
					int rj = order[i];
					if(reads[rj].chrIdx != ci || reads[rj].beg != pos) break;
					sum += parseWeight(reads, rj);
					i++;
				}
				printf("%s\t%d\t%d\t.\t%g\t+\n", chroms[ci].c_str(), pos, pos+1, sum);
			}
		}
		else
		{
			const bool fullLine = fRal || useMmap;
			for(int i = 0; i < totalReads; i++)
			{
				int ri = order[i];
				printRead(reads, ri, chroms[reads[ri].chrIdx].c_str(), fullLine, stdout);
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
			int* chromTable = (int*) calloc((size_t)maxChrLen + 1, sizeof(int));
			int* numReadsBeg = (sortMode == 'b')
				? (int*) calloc((size_t)maxChrLen + 1, sizeof(int))
				: NULL;
			for(const std::string& chrom : chroms)
			{
				cerr << "Sorting " << chrom << endl;
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
			cerr << "Parallel bucket sort using up to " << numThreads << " threads ("
			     << chroms.size() << " chromosomes)" << endl;
			std::vector<int> idxs(chroms.size());
			for(size_t i = 0; i < chroms.size(); i++) idxs[i] = (int)i;

			std::mutex printMtx;
			std::condition_variable printCv;
			int nextChromToPrint = 0;

			std::for_each(std::execution::par, idxs.begin(), idxs.end(),
				[&](int ci) {
					const std::string& chrom = chroms[ci];
					string2chrInfoT::iterator cit = chrInfo.find(chrom);
					int thisChrLen = cit->second.len;

					char* memBuf = NULL;
					size_t memBufSz = 0;
					FILE* sink = open_memstream(&memBuf, &memBufSz);
					if(!sink)
					{
						cerr << "open_memstream failed for chrom " << chrom << endl;
						exit(1);
					}
					{
						std::vector<int> chromTable((size_t)thisChrLen + 1, 0);
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
	time(&tend);
	double difftime_sorting = difftime(tend, tstart);
	fprintf(stderr, "Sorting has taken %d seconds\n", (int)(difftime_sorting+0.5));
	return 0;
}

/* some stats:
 * for 108'170'237 reads required 20GB Reading has taken 241 seconds
 * Sorting has taken 834 seconds
 *
 * After optimizing the same set:
 * Reading has taken 113 seconds
 * Sorting has taken 30 seconds
 * 9.7 GB RAM
 */
