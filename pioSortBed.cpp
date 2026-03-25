#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <parallel/algorithm>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <unordered_map>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "CLI11.hpp"

/* Copyright: Piotr Balwierz */

using namespace std;

// Field layout minimizes padding: pointer first, then ints, then char.
// 24 bytes per read on 64-bit (8 + 4 + 4 + 4 + 1 + 3 padding).
// Removing str would NOT shrink the struct (still 24 bytes due to pointer alignment).
//
// The `next` field is used by bucket sort to thread per-chromosome linked lists.
// The `chrIdx` field is used by the classic sort path to store a small integer
// chromosome index for fast comparison.  They occupy the same 4 bytes (union)
// because the two sort paths never need both simultaneously:
//   - During parsing, `next` is populated (linked lists).
//   - If classic sort is chosen, we walk each list once to stamp `chrIdx`,
//     overwriting `next` which is no longer needed.
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

// Possible issues: UTF8 symbols, long lines;
const int lineBufSize = 1024;		// buffer for the read file lines (stdin path only; mmap has no line length limit)
const int chrLenLimit = 1000000000;  // 500Mb. we don't believe there are longer chromosomes;
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
// Used for stdin input (all lines) and for collapse mode (weight strings).
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
	size_t      chunkSize;  // size to use when adding a new chunk

	Arena(size_t initialCap, size_t chunkSz) : chunkSize(chunkSz)
	{
		head = new ArenaChunk(initialCap);
	}
	~Arena()
	{
		ArenaChunk* c = head;
		while (c) { ArenaChunk* nx = c->next; delete c; c = nx; }
	}

	// Copy len bytes from src into the arena and return a pointer to the copy.
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
// Skip over a field (tab or space delimited), advancing *pp past the delimiter.
static inline bool skipField(const char** pp)
{
	const char* p = *pp;
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	if (*p == '\t' || *p == ' ') { *pp = p + 1; return true; }
	*pp = p;
	return false;
}

// Parse a non-negative integer, advancing *pp past the digits.
static inline int parseUInt(const char** pp)
{
	const char* p = *pp;
	if (*p < '0' || *p > '9') return -1;
	int val = 0;
	while (*p >= '0' && *p <= '9') { val = val * 10 + (*p - '0'); p++; }
	*pp = p;
	return val;
}

// Copy field text into dst (up to dstMax-1 chars), NUL-terminate, advance *pp.
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

// Minimal BED parser: only chr (as pointer+length into buf), beg, end.
// Also returns tailPtr: pointer to everything after end field (fields 4+).
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

// Full BED parser: chr (as pointer+length), beg, end, weight, strand.
// Also returns tailPtr (same as parseBedLine3).
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
	*tailPtr = p;  // capture tail before parsing deeper fields

	if (*p != '\t' && *p != ' ') return 3;
	p++;

	// field 3: id — skip it
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	if (*p != '\t' && *p != ' ') return 3;
	p++;

	// field 4: weight
	int n = copyField(&p, weightBuf, weightMax);
	if (n == 0) return 3;

	// field 5: strand
	if (*p != '\t' && *p != ' ') return 4;
	p++;
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

float sumList(seqread* reads, int* readBuff, int foo);
float sumList2(seqread*, int*, int);

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

// Low-memory file mode (SSD-friendly):
// Pass 1: walk mmap once, build per-chromosome linked lists of line offsets.
// Pass 2: process one chromosome at a time (parse, sort, print), so peak RAM
// depends on the largest chromosome chunk, not the whole file.
static int lowMemSortMmap(char* mmapBase, size_t mmapSize,
                          bool fRal, int fCollapse, char sortMode, int numThreads)
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
	std::sort(chroms.begin(), chroms.end());

	time(&tend);
	fprintf(stderr, "Reading has taken %d seconds\n", (int)(difftime(tend, tstart)+0.5));
	cerr << "We have " << nodeCount - 1 << " regions.\n"
		 << chroms.size() << " chromosomes\nSorting..." << endl;
	time(&tstart);

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
				if(!sscanf(weight, "%f", &r.weight))
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
			__gnu_parallel::sort(recs.begin(), recs.end(), cmp);

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

// Parsing loop, templated on UseMmap to eliminate per-line branch.
// Uses parseBedLine3 (3 fields only) for the common BED case,
// parseBedLineFull when weight/strand are needed, parseRalLine for RAL.	// For BED, stores only the "tail" (fields 4+) in reads[].line, saving ~50% memory.
	// For RAL, stores the full line (RAL output format differs from BED).// Same-chr check uses length-first comparison to avoid string construction.
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
			// Strip trailing \r\n now so parsers and tailPtr see clean data
			size_t len = strlen(linePtr);
			if (len >= 2 && linePtr[len-2] == '\r' && linePtr[len-1] == '\n')
				{ linePtr[len-2] = '\0'; }
			else if (len >= 1 && linePtr[len-1] == '\n')
				{ linePtr[len-1] = '\0'; }
		}

		if(readCount == currMaxReads)
		{
			currMaxReads = (int) (currMaxReads * 2);
			reads = (seqread*) realloc(reads, currMaxReads * sizeof(seqread));
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
				// RAL: store full line (output format differs from BED)
				if(UseMmap)
					reads[readCount].line = linePtr;
				else
					reads[readCount].line = arena->alloc(linePtr, strlen(linePtr) + 1);
			}
			else if(!fCollapse)
			{
				if(UseMmap)
				{
					// mmap: point to full line (zero-copy, fast output)
					reads[readCount].line = linePtr;
				}
				else
				{
					// stdin: store only tail to save arena memory
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
			// Avoids constructing a std::string for the common case (same chr).
			if((int)thisChrIt->first.size() != chrLen ||
			   memcmp(thisChrIt->first.data(), chrPtr, chrLen) != 0)
			{	// different chromosome
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
			{	// same chromosome as previous line
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

int main(int argc, char *argv[])
{
	// 1MB output buffer
	char* outBuf = (char*) malloc(1 << 20);
	setvbuf(stdout, outBuf, _IOFBF, 1 << 20);

	CLI::App app{"Ultra fast bed file sorter\nPiotr Balwierz, 2012\n\n"
		"Results are equivalent to \"LC_ALL=C sort -k1,1 -k2,2n file.bed\"\n"
		"or \"sort -k1,1 -k2,2n -k3,3n file.bed\" if --sort=b enabled.\n"
		"Uses one thread by default and sorts at ~disk IO throughput limits.\n\n"
		"Input file should contain a new line character in the end of the last line.\n\n"
		"Compilation-time limits:\n"
		"  Line length limit: " + to_string(lineBufSize) + " (stdin only; no limit for file input)\n"
		"  Chromosome name limit: " + to_string(chrNameBufSize) + "\n"
		"  Chromosome length limit: " + to_string(chrLenLimit/1000000) + "Mbp\n"
		"  Chromosome/contig number limit: unlimited\n"
		"  Read number limit: unlimited (but will be kept in the memory)\n"};
	app.set_version_flag("-V,--version", "2.0.0");

	string inputFile;
	char sortMode = 's';
	int fCollapse = 0;
	int fRal = 0;
	int lowMemSSD = 0;
	int bucketCutoff = defaultBucketCutoff;
	int numThreads = 0;

	app.add_option("input-file", inputFile, "input file; put \"-\" to read from stdin")
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

	CLI11_PARSE(app, argc, argv);

	if(numThreads <= 0)
		numThreads = omp_get_max_threads();
	omp_set_num_threads(numThreads);

	int readCount = 1;	// count == 0 is used as an end of a list.
	int currMaxReads = 1024;
	long fileSize = -1; // -1 means stdin / unknown

	// I/O state: either mmap (file) or FILE* (stdin)
	char* mmapBase = NULL;
	size_t mmapSize = 0;
	bool useMmap = false;
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
			close(fd);  // fd can be closed after mmap
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

			// Arena only needed in collapse mode (for weight strings)
			if(fCollapse)
			{
				size_t arenaInit = (size_t)(fileSize / 5);
				if(arenaInit < 4096) arenaInit = 4096;
				arena = new Arena(arenaInit, 64UL * 1024 * 1024);
			}
		}

		if(lowMemSSD)
		{
			if(!useMmap)
			{
				cerr << "Error: --low-mem-ssd requires file input (not stdin)" << endl;
				return 1;
			}
			return lowMemSortMmap(mmapBase, mmapSize, fRal, fCollapse, sortMode, numThreads);
		}

		// allocate the guessed amount of memory:
		seqread *reads = (seqread*) malloc(currMaxReads * sizeof(seqread));

		time_t tstart, tend;
		time(&tstart);

		int maxChrLen = 0;
		string2chrInfoT chrInfo;
		pair<string2chrInfoT::iterator,bool> insResult =
			chrInfo.insert(make_pair(string(""), chrInfoT()));
		string2chrInfoT::iterator thisChrIt = insResult.first;

		// Parsing loop (templated on mmap vs stdin to eliminate per-line branch)
		if(useMmap)
			parseLines<true>(mmapBase, mmapSize, fh,
			                 reads, readCount, currMaxReads,
			                 chrInfo, thisChrIt, arena,
			                 fRal, fCollapse, sortMode);
		else
			parseLines<false>(mmapBase, mmapSize, fh,
			                  reads, readCount, currMaxReads,
			                  chrInfo, thisChrIt, arena,
			                  fRal, fCollapse, sortMode);
		if(fh) fclose(fh);
		time(&tend);
		double difftime_reading = difftime(tend, tstart);
		fprintf(stderr, "Reading has taken %d seconds\n", (int)(difftime_reading+0.5));

		/***********************************
		 *  the actual sorting starts here *
		 * *********************************/
		chrInfo.erase("");
		// list all the chromosome NAMES here:
		int nChrom = chrInfo.size();
		std::vector<std::string> chroms;
		for(string2chrInfoT::iterator it=chrInfo.begin(); it!=chrInfo.end(); it++)
		{
			cerr << it->first << ": " << it->second.len << endl;
			chroms.push_back(it->first);
			if(maxChrLen < it->second.len)
			{
				maxChrLen = it->second.len;
			}
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
		std::sort(chroms.begin(), chroms.end());

		// Assign sequential chromosome indices (alphabetical order) so the
		// classic sort path can compare chromosomes with a single int compare
		// instead of a string compare.
		for(int ci = 0; ci < (int)chroms.size(); ci++)
			chrInfo.find(chroms[ci])->second.idx = ci;

		// Choose sort strategy: classic O(n log n) for small files,
		// bucket/counting sort O(n + m) for large files.
		// The bucket sort allocates chromTable[maxChrLen+1] which can be up to 4 GB.
		// For small files this is wasteful — the classic path uses only
		// readCount * 4 bytes for the index array, and std::sort with an inlined
		// comparator is very efficient for moderate n.
		int totalReads = readCount - 1;
		// bucketCutoff == 0 means "always use bucket sort" (skip the classic path).
		// Otherwise use bucket sort when we have at least bucketCutoff reads.
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
			 * We sort an index array rather than the 24-byte seqread structs.
			 * Swapping 4-byte indices is cheaper and keeps the reads array
			 * in place (important when line pointers reference the mmap buffer).
			 *
			 * Steps:
			 *  1. Walk each chromosome's linked list to stamp chrIdx and
			 *     collect all read indices into a flat array.
			 *  2. std::sort the index array with an inlined lambda comparator.
			 *  3. Linear scan to print (and optionally collapse).
			 *************************************************************/

			// Step 1: Walk linked lists, stamp chrIdx, build index array.
			// The `next` field is no longer needed after this — we reuse its
			// storage for chrIdx (they share a union).
			int* order = (int*) malloc(totalReads * sizeof(int));
			int oi = 0;
			for(int ci = 0; ci < (int)chroms.size(); ci++)
			{
				string2chrInfoT::iterator cit = chrInfo.find(chroms[ci]);
				int curr = cit->second.lastRead;
				while(curr)
				{
					int nxt = reads[curr].next;  // save before overwriting
					reads[curr].chrIdx = ci;
					order[oi++] = curr;
					curr = nxt;
				}
			}

			// Step 2: Sort the index array.
			// When numThreads == 1, use std::sort (introsort) which inlines the
			// comparator lambda — best single-threaded performance.
			// When numThreads > 1, use __gnu_parallel::sort which splits the
			// work across OpenMP threads.  The parallel version cannot inline
			// the comparator as aggressively, so for small inputs the overhead
			// of thread management may negate the parallelism gains.  For
			// medium-to-large inputs (1M+ reads) the speedup is significant.
			//
			// Primary key: chromosome index (int compare, not string).
			// Secondary key: sort position (beg, or 5' end).
			// Tertiary key (--sort b only): end position.
			if(numThreads == 1)
			{
				// Single-threaded: std::sort with inlined lambda — fastest for 1 core.
				if(sortMode == 'b')
				{
					std::sort(order, order + totalReads, [reads](int a, int b) {
						if(reads[a].chrIdx != reads[b].chrIdx) return reads[a].chrIdx < reads[b].chrIdx;
						if(reads[a].beg != reads[b].beg) return reads[a].beg < reads[b].beg;
						return reads[a].end < reads[b].end;
					});
				}
				else if(sortMode == '5')
				{
					std::sort(order, order + totalReads, [reads](int a, int b) {
						if(reads[a].chrIdx != reads[b].chrIdx) return reads[a].chrIdx < reads[b].chrIdx;
						int posA = (reads[a].str == '-') ? reads[a].end : reads[a].beg;
						int posB = (reads[b].str == '-') ? reads[b].end : reads[b].beg;
						return posA < posB;
					});
				}
				else
				{
					std::sort(order, order + totalReads, [reads](int a, int b) {
						if(reads[a].chrIdx != reads[b].chrIdx) return reads[a].chrIdx < reads[b].chrIdx;
						return reads[a].beg < reads[b].beg;
					});
				}
			}
			else
			{
				// Multi-threaded: __gnu_parallel::sort partitions the array
				// across OpenMP threads.  Each partition is sorted independently,
				// then merged — giving ~N× speedup on the sort phase alone.
				cerr << "Parallel sort using " << numThreads << " threads" << endl;
				if(sortMode == 'b')
				{
					__gnu_parallel::sort(order, order + totalReads, [reads](int a, int b) {
						if(reads[a].chrIdx != reads[b].chrIdx) return reads[a].chrIdx < reads[b].chrIdx;
						if(reads[a].beg != reads[b].beg) return reads[a].beg < reads[b].beg;
						return reads[a].end < reads[b].end;
					});
				}
				else if(sortMode == '5')
				{
					__gnu_parallel::sort(order, order + totalReads, [reads](int a, int b) {
						if(reads[a].chrIdx != reads[b].chrIdx) return reads[a].chrIdx < reads[b].chrIdx;
						int posA = (reads[a].str == '-') ? reads[a].end : reads[a].beg;
						int posB = (reads[b].str == '-') ? reads[b].end : reads[b].beg;
						return posA < posB;
					});
				}
				else
				{
					__gnu_parallel::sort(order, order + totalReads, [reads](int a, int b) {
						if(reads[a].chrIdx != reads[b].chrIdx) return reads[a].chrIdx < reads[b].chrIdx;
						return reads[a].beg < reads[b].beg;
					});
				}
			}

			// Step 3: Print in sorted order.
			if(fCollapse)
			{
				// Collapse mode: adjacent reads at the same (chr, position) are
				// merged by summing weights.  After sorting, identical positions
				// are contiguous, so a single linear scan suffices.
				int i = 0;
				while(i < totalReads)
				{
					int ri = order[i];
					int ci = reads[ri].chrIdx;
					int pos = reads[ri].beg;
					float sum = 0.0;
					// Accumulate all reads at the same (chr, pos)
					while(i < totalReads)
					{
						int rj = order[i];
						if(reads[rj].chrIdx != ci || reads[rj].beg != pos) break;
						float w;
						if(!sscanf(reads[rj].line, "%f", &w))
						{
							cerr << "Malformed weight " << reads[rj].line << endl;
							exit(1);
						}
						sum += w;
						i++;
					}
					printf("%s\t%d\t%d\t.\t%g\t+\n", chroms[ci].c_str(), pos, pos+1, sum);
				}
			}
			else
			{
				if(fRal || useMmap)
				{
					// Full line pointer (RAL or mmap): fast single fputs
					for(int i = 0; i < totalReads; i++)
					{
						fputs_unlocked(reads[order[i]].line, stdout);
						fputc_unlocked('\n', stdout);
					}
				}
				else
				{
					// stdin BED: reconstruct from tail
					for(int i = 0; i < totalReads; i++)
					{
						int ri = order[i];
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
			 * Allocates chromTable[maxChrLen+1] and scatters reads by position.
			 * Wins for large files where n >> log(n), making the linear scan
			 * over chromosome length worthwhile.
			 *************************************************************/
			int* chromTable = (int*)  calloc(maxChrLen+1, sizeof(int));
			// numReadsBeg only needed for --sort b mode
			int* numReadsBeg = NULL;
			if(sortMode == 'b')
				numReadsBeg = (int*) calloc(maxChrLen+1, sizeof(int));

			for(std::vector<std::string>::iterator it=chroms.begin(); it!=chroms.end(); it++)
			{
				cerr << "Sorting " << *it << endl;
				/** PART 1 of sorting
				 * building an array of lists
				 **/
				// Single hash lookup per chromosome
				string2chrInfoT::iterator cit = chrInfo.find(*it);
				int thisChrLen = cit->second.len;
				int currRead = cit->second.lastRead;
				int maxNumReads = 0;	// the maximum number of reads at any position
				while(currRead)
				{
					int chosenPosition = (sortMode == '5') ? (reads[currRead].str == '-' ? reads[currRead].end : reads[currRead].beg) : reads[currRead].beg;
					int oldPtr = chromTable[chosenPosition];
					chromTable[chosenPosition] = currRead;
					currRead = reads[currRead].next;		// this and the next line cannot be exchanged
					reads[chromTable[chosenPosition]].next = oldPtr;
					if(sortMode == 'b')
					{
						numReadsBeg[chosenPosition] ++;
						if(numReadsBeg[chosenPosition] > maxNumReads)
							maxNumReads ++;	// it can only increase by one.
					}
				}
				/** PART 2 of sorting:
				 * printing and destroying the count list:
				 **/
				if(sortMode == 'b')
				{
					// we want to sort by end position too.
					int* readBuff = (int*) calloc(maxNumReads, sizeof(int));
					for(int pos=0; pos<=thisChrLen; pos++)
					{
						if(chromTable[pos])
						{
							int i=0;
							int foo = numReadsBeg[pos];
							while(chromTable[pos])
							{
								readBuff[i++] = chromTable[pos];
								chromTable[pos] = reads[chromTable[pos]].next;
							}
							// std::sort with lambda inlines the comparator (faster than qsort function pointer)
							std::sort(readBuff, readBuff + foo, [reads](int a, int b) {
								return reads[a].end < reads[b].end;
							});
							// print it out:
							if(fCollapse)
								printf("%s\t%d\t%d\t.\t%g\t+\n", it->c_str(), pos, pos+1, sumList(reads, readBuff, foo));
							else if(fRal || useMmap)
								for(int j = 0; j < foo; ++j)
								{
									fputs_unlocked(reads[readBuff[j]].line, stdout);
									fputc_unlocked('\n', stdout);
								}
							else
								for(int j = 0; j < foo; ++j)
									writeBedLine(it->c_str(), reads[readBuff[j]].beg,
									             reads[readBuff[j]].end,
									             reads[readBuff[j]].line, stdout);
							numReadsBeg[pos] = 0;
						}
					}
					free(readBuff);
				}
				else  // sortMode in {5, s}
				{
					// we don't want to sort by the end position: simply print what is attached to any position in the genome
					for(int pos=0; pos<=thisChrLen; pos++)
						if(chromTable[pos])
						{
							if(fCollapse)
								printf("%s\t%d\t%d\t.\t%g\t+\n", it->c_str(), pos, pos+1, sumList2(reads, chromTable, pos));
							else if(fRal || useMmap)
								while(chromTable[pos])
								{
									fputs_unlocked(reads[chromTable[pos]].line, stdout);
									fputc_unlocked('\n', stdout);
									chromTable[pos] = reads[chromTable[pos]].next;
								}
							else
								while(chromTable[pos])
								{
									int ri = chromTable[pos];
									writeBedLine(it->c_str(), reads[ri].beg, reads[ri].end,
									             reads[ri].line, stdout);
									chromTable[pos] = reads[ri].next;
								}
						}
				}
			} // end per-chromosome loop (bucket sort)
			free(chromTable);
			if(numReadsBeg) free(numReadsBeg);
		} // end bucket sort path
		time(&tend);
		double difftime_sorting = difftime(tend, tstart);
		fprintf(stderr, "Sorting has taken %d seconds\n", (int)(difftime_sorting+0.5));
		return 0;
}

float sumList(seqread* reads, int* readBuff, int foo)
{
	float sum = 0.0;
	for(int j = 0; j < foo; ++j)
	{
		float w;
		if(! sscanf(reads[readBuff[j]].line, "%f", &w))
		{
			cerr << "Malformed weight " << reads[readBuff[j]].line << endl;
			exit(1);
		}
		sum += w;
	}
	return sum;
}

float sumList2(seqread* reads, int* chromTable, int pos)
{
	float sum = 0.0;
	while(chromTable[pos])
	{
		float w;
		if(! sscanf(reads[chromTable[pos]].line, "%f", &w))
		{
			cerr << "Malformed weight " << reads[chromTable[pos]].line << endl;
			exit(1);
		}
		sum += w;
		chromTable[pos] = reads[chromTable[pos]].next;
	}
	return sum;
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
