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
const int defaultBucketCutoff = 10000000;


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

// BED line parser.
static int parseBedLine(const char* buf,
                        char* chr, int chrMax,
                        int* beg, int* end,
                        char* weightBuf, int weightMax,
                        char* strandChar)
{
	const char* p = buf;

	int n = copyField(&p, chr, chrMax);
	if (n == 0) return 0;

	int v = parseUInt(&p);
	if (v < 0) return 1;
	*beg = v;
	if (*p == '\t' || *p == ' ') p++;

	v = parseUInt(&p);
	if (v < 0) return 2;
	*end = v;

	if (*p != '\t' && *p != ' ') return 3;
	p++;

	// field 3: id — skip it
	while (*p && *p != '\t' && *p != ' ' && *p != '\r' && *p != '\n') p++;
	if (*p != '\t' && *p != ' ') return 3;
	p++;

	// field 4: weight
	n = copyField(&p, weightBuf, weightMax);
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

	string inputFile;
	char sortMode = 's';
	int fCollapse = 0;
	int fRal = 0;
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

		// allocate the guessed amount of memory:
		seqread *reads = (seqread*) malloc(currMaxReads * sizeof(seqread));

		time_t tstart, tend;
		time(&tstart);

		int maxChrLen = 0;
		string2chrInfoT chrInfo;
		pair<string2chrInfoT::iterator,bool> insResult =
			chrInfo.insert(make_pair(string(""), chrInfoT()));
		string2chrInfoT::iterator thisChrIt = insResult.first;

		// Parsing loop
		char* mmapCur = mmapBase;
		char* mmapLim = mmapBase + mmapSize;
		char stdinBuf[lineBufSize];

		while(1)
		{
			char* linePtr;

			if(useMmap)
			{
				if(mmapCur >= mmapLim) break;
				linePtr = mmapCur;
				// memchr uses SIMD on modern glibc — much faster than fgets char-by-char copy
				char* nl = (char*) memchr(mmapCur, '\n', (size_t)(mmapLim - mmapCur));
				if(nl)
				{
					if(nl > mmapCur && *(nl - 1) == '\r') *(nl - 1) = '\0';
					*nl = '\0';
					mmapCur = nl + 1;
				}
				else
				{
					// Last line without trailing newline.
					// The extra byte mapped past EOF is guaranteed to be zero (NUL sentinel).
					mmapCur = mmapLim;
				}
			}
			else
			{
				if(!fgets_unlocked(stdinBuf, lineBufSize, fh)) break;
				linePtr = stdinBuf;
			}

			if(readCount == currMaxReads)
			{
				currMaxReads = (int) (currMaxReads * 2);
				reads = (seqread*) realloc(reads, currMaxReads * sizeof(seqread));
			}

			int beg = 0;
			int end = 0;
			char chr[chrNameBufSize];
			char strandChar = '+';
			char weight[chrNameBufSize];
			weight[0] = '0'; weight[1] = '\0'; // default weight

			int numArgsRead;
			if (fRal)
				numArgsRead = parseRalLine(linePtr, chr, chrNameBufSize,
				                           &beg, &end,
				                           weight, chrNameBufSize,
				                           &strandChar);
			else
				numArgsRead = parseBedLine(linePtr, chr, chrNameBufSize,
				                           &beg, &end,
				                           weight, chrNameBufSize,
				                           &strandChar);

			if(numArgsRead >= 3)
			{
				if(beg < 0)
				{
					cerr << "Error: negative coordinates in the bed file\n" << linePtr << endl;
					exit(1);
				}
				if(!fCollapse)
				{
					if(useMmap)
					{
						// Line is already NUL-terminated in the mmap buffer (no copy needed)
						reads[readCount].line = linePtr;
					}
					else
					{
						// Strip trailing \r\n or \n before storing
						size_t len = strlen(linePtr);
						if (len >= 2 && linePtr[len-2] == '\r' && linePtr[len-1] == '\n')
						{
							linePtr[len-2] = '\0';
							len -= 2;
						}
						else if (len >= 1 && linePtr[len-1] == '\n')
						{
							linePtr[len-1] = '\0';
							len--;
						}
						reads[readCount].line = arena->alloc(linePtr, len + 1);
					}
				}
				else
				{
					// In collapse mode store just the weight string.
					size_t wlen = strlen(weight);
					reads[readCount].line = arena->alloc(weight, wlen + 1);
				}
				reads[readCount].beg = beg;
				reads[readCount].end = end;
				reads[readCount].str = strandChar;
				// depending on parameters choose the right position:
				int chosenPosition = (sortMode == '5') ? (strandChar == '-' ? end : beg) : beg;

				// no lookup if we deal with the same chromosome as in the previous line:
				if(thisChrIt->first.compare(chr))
				{	// it is a different chromosome :-(
					pair<string2chrInfoT::iterator,bool> ins =
						chrInfo.insert(make_pair(string(chr), chrInfoT()));
					thisChrIt = ins.first;
					if(!ins.second)
					{
						// We already have seen reads from this chromosome:
						if(thisChrIt->second.len < chosenPosition) thisChrIt->second.len = chosenPosition;
						// push head to the list:
						reads[readCount].next = thisChrIt->second.lastRead;
					}
					else
					{
						// This chromosome is a new one
						// create the end element in the list:
						reads[readCount].next = 0;
						thisChrIt->second.len = chosenPosition;
					}
				}
				else
				{	// the same chromosome as in the previous line
					if(thisChrIt->second.len < chosenPosition) thisChrIt->second.len = chosenPosition;
					// push head to the list:
					reads[readCount].next = thisChrIt->second.lastRead;
				}
				thisChrIt->second.lastRead = readCount;
			}
			else
			{
				cerr << "Error in parsing line: " << linePtr << endl
					<< "Perhaps this line is malformed?" << endl;
				return(1);
			}
			readCount ++;
		} // while over the lines in the file
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
				for(int i = 0; i < totalReads; i++)
				{
					fputs_unlocked(reads[order[i]].line, stdout);
					fputc_unlocked('\n', stdout);
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
							else
								for(int j = 0; j < foo; ++j)
								{
									fputs_unlocked(reads[readBuff[j]].line, stdout);
									fputc_unlocked('\n', stdout);
								}
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
							else
								while(chromTable[pos])
								{
									fputs_unlocked(reads[chromTable[pos]].line, stdout);
									fputc_unlocked('\n', stdout);
									chromTable[pos] = reads[chromTable[pos]].next;
									// no deallocation here for speed
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
