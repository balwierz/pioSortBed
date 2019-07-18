#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

/* Copyright: Piotr Balwierz */

/* Correct way of compiling it is:
 * /import/bc2/soft/app/gcc/4.6.3/Linux/bin/g++ -I/import/bc2/soft/app/boost/1.42.0/Linux/include -L/import/bc2/soft/app/boost/1.42.0/Linux/lib pioSortBed6.cpp -o pioSortBed -O3 -lboost_program_options -static
 * */



using namespace std;

// we need 32 bit to store the chromosome coordinate. On 64bit opteron linux int is 4bytes, int_fast32_t is however 8bytes
// and we don't want to store such numbers;
class seqread
{
	public:
	int beg;
	int end;
	char str;
	char* line; // keeps the full line; in case of collapse option just the weight as string.
	int next;
};

// make sure we call find or [] on a chr only once!
class chrInfoT
{
	public:
	int len;		// length (the max beg position)
	int lastRead;	// keeps the head of the list of all reads at a given chromosome
	chrInfoT() : len(0), lastRead(0) {}
};

typedef boost::unordered_map<string, chrInfoT> string2chrInfoT;

// Possible issues: UTF8 symbols, long lines;
const int lineBufSize = 1024;		// buffer for the read file lines
const int chrLenLimit = 1000000000;  // 500Mb. we don't believe there are longer chromosomes;
									// Just a safety check, not to allocate TBytes of ram
									// Increase if you know what you are doing.
									// And change the text chrTooLongMsg below.
									// Corresponds to 1000M * 4B = 4GB on a 64bit machine
const char chrTooLongMsg[] = "That is more than 3 times the length of human chr1! If you are sure this is correct recompile the program with an increased chrLenLimit value";
const int chrNameBufSize = 256;		// buffer for chromosme names like "chr1_random"
									// and strand too, although it should be only "+" or "-"


int estimateReadNumber(FILE*);
int cmpString(const void*, const void*);
int cmpSeqreadPtrEnd(const void* a, const void* b);
int cmpReadIndxEnd(const void* a, const void* b);
float sumList(seqread* reads, int* readBuff, int foo);
float sumList2(seqread*, int*, int);
using namespace boost;
seqread* readsGlobal; // a hack for sorting. Unfortunately we need a global

namespace po = boost::program_options;
int main(int argc, char *argv[])
{
	//cout << sizeof(int) << " " << sizeof(int32_t) << " " << sizeof(int_fast32_t) << endl; return 0;
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
		    ("help,h", "produce this help message")
		    ("ral,r", "input is in RAL format")
		    ("sort,s", po::value<char>()->default_value('s'), "arg can be: s-start, b-both, 5-5'end\n"
						"s: start coordinate [default]\n"
						"b: sort also by the end coordinate, aka \"LC_ALL=C sort -k1,1 -k2,2n -k3,3n\". Don't use if you don't need it, as it requires more computation.\n"
						"5: sort by 5' -p- bond (2nd column for + strand and 3rd column for - strand)\n"
			)
		    ("collapse", "collapse BEDWEIGHT by summing weights of the reads and truncate the coordinates, set strand to \"+,\" id to \".\". Makes no sense for --sort=b")
		    ("input-file", po::value<string>(), "input file; put \"-\" to read from standard input (e.g. cat file.bed | pioSortBed)")
		;
		
		po::positional_options_description p;
		p.add("input-file", -1);
		
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
	
		if (vm.count("help") || ! vm.count("input-file"))
		{
		    cerr << "Ultra fast bed file sorter\nPiotr Balwierz, 2012\n\n"
			<< "Usage: piotrSortBed file.bed > outfile.bed\n\n"
			<< "Results are equivalent to \"LC_ALL=C sort -k1,1 -k2,2n file.bed\"\n"
			<< "or \"sort -k1,1 -k2,2n -k3,3n file.bed\" if --end option enabled.\n"
			<< "Uses one thread and sorts at ~disk IO throughput limits.\n\n"
			<< "Input file should contain a new line character in the end of the last line.\n\n"
			<< "Compilation-time limits:\nLine length limit: " << lineBufSize
			<< "\nChromosome name limit: " << chrNameBufSize
			<< "\nChromosome length limit: " << chrLenLimit/1000000 << "Mbp"
			<< "\nChromosome/contig number limit: unlimited"
			<< "\nRead number limit: unlimited (but will be kept in the memory)\n" << endl
		    << desc << "\n";
		    return 1;
		}
		//if (vm.count("input-file"))
		//{
		//	cout << "Input file: " << vm["input-file"].as<string>() << "\n";
		//}
		char* lineFormat;
		char lineFormatBed[] = "%s %d %d %*s %s %s\n";
		char lineFormatRal[] = "%*s %s %*s %d %d %s %*s %s\n";  //sq2  chr2L   +  12913266   12913300    0   TTCC  ATCC
		if(vm.count("ral"))
		{	
			lineFormat = lineFormatRal;
		}
		else
		{
			lineFormat = lineFormatBed;
		}
			
	
		int readCount = 1;	// count == 0 is used as an end of a list.
		FILE* fh = NULL;
		int currMaxReads = 1024;
		if(vm["input-file"].as<string>() == "-")
		{
			fh = stdin;
			cerr << "Reading data from standard input" << endl;
		}
		else
		{
			fh = fopen(vm["input-file"].as<string>().c_str(), "r");
			if(fh == NULL)
			{
				cerr << "Error opening " << vm["input-file"].as<string>() << endl;
				return 1;
			}
			// allocate hopefully enough of memory for the storage of reads:
			// we estimate the number of reads in the file and add a 10% margin
			currMaxReads = (int)(estimateReadNumber(fh) * 1.1);
		}
		// these are to not do lookups in the hash each time:
		int fCollapse = vm.count("collapse");
		char sortMode = vm["sort"].as<char>();
		
		// allocate the guessed amount o memory:
		seqread *reads = (seqread*) malloc(currMaxReads * sizeof(seqread));
		
		time_t start, end;
		time(&start);
		
		int maxChrLen = 0;
		string2chrInfoT chrInfo;
		chrInfo[""]; // create a dummy chromosome
		string2chrInfoT::iterator thisChrIt = chrInfo.find("");
		while(1)
		{
			char buff[lineBufSize];
			if(fgets(buff, lineBufSize, fh))
			{
				if(readCount == currMaxReads)
				{
					// expand the size by 30% 
					currMaxReads = (int) (currMaxReads * 1.3);
					reads = (seqread*) realloc(reads, currMaxReads * sizeof(seqread));
					// cerr << "Reallocating reads\n";
				}
				// parse the chrom, strand and begining (but not more)
				int beg;
				int end;
				char chr[chrNameBufSize];
				char str[chrNameBufSize];	// in ral str means simply the last colum == genomic seq
				char weight[chrNameBufSize];
				int numArgsRead = sscanf(buff, lineFormat, chr, &beg, &end, weight, str);
				
				if(numArgsRead >= 3)
				{
					// we have read to the end of line
					if(beg < 0)
					{
						cerr << "Error: negative coordinates in the bed file\n" << buff << endl;
						exit(1);
					}
					if(!fCollapse)
					{
						reads[readCount].line = (char*) malloc(strlen(buff)+1);
						strcpy(reads[readCount].line, buff);
					}
					else
					{
						reads[readCount].line = (char*) malloc(strlen(weight)+1);
						strcpy(reads[readCount].line, weight);
					}
					reads[readCount].beg = beg;
					reads[readCount].end = end;
					reads[readCount].str = str[0];
					// depending on parameters choose the right position:
					int chosenPosition = (sortMode == '5') ? (strcmp(str, "+") ? end : beg) : beg;
	
					// no lookup if we deal with the same chromosome as in the previous line:
					if(thisChrIt->first.compare(chr))
					{	// it is a different chromosome :-(	
						thisChrIt = chrInfo.find(chr);
						if(thisChrIt != chrInfo.end())
						{
							// We alredy have seen read from this chromosome: 
							if(thisChrIt->second.len < chosenPosition) thisChrIt->second.len = chosenPosition;
							// push head to the list:
							reads[readCount].next = thisChrIt->second.lastRead;
						}
						else
						{
							// This chromosome is a new one
							//cerr << "Initializing " << chr << endl;
							// create the end element in the list:
							reads[readCount].next = 0;
							// the iterator is invalidated now; fix it:
							chrInfo[chr];
							thisChrIt = chrInfo.find(chr);
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
					cerr << "Error in parsing line: " << buff << endl
						<< "Perhaps this line is malformed or too long "
						<< "(not enough buffer of " << lineBufSize << " bytes)?" << endl;
					return(1);
				}
			}
			else
			{
				break; // the last line was read
			}
			readCount ++;
		} // while over the lines in the file
		fclose(fh);
		time(&end);
		double difftime_reading = difftime(end, start);
		fprintf(stderr, "Reading has taken %d seconds\n", (int)(difftime_reading+0.5));
		
		/***********************************
		 *  the actual sorting starts here *
		 * *********************************/
		cerr << "We have " << readCount-1 << " regions. Sorting..." << endl;
		chrInfo.erase("");
		// list all the chromosome NAMES here:
		int nChrom = chrInfo.size();
		string* chroms = new string[nChrom];
		int chrI = 0;
		for(string2chrInfoT::iterator it=chrInfo.begin(); it!=chrInfo.end(); it++)
		{
			chroms[chrI++] = it->first;
			if(maxChrLen < it->second.len)
			{
				maxChrLen = it->second.len;
			}
		}
		//cout << "Longest chr is " << maxChrLen << endl;
		if(maxChrLen > chrLenLimit)
		{
			cerr << "Error: There is a region starting at " << maxChrLen << "." << endl
				<<  chrTooLongMsg << endl;
			exit(1);
		}
		time(&start);
		qsort(chroms, nChrom, sizeof(string), cmpString);
		int* chromTable = (int*)  calloc(maxChrLen+1, sizeof(int)); // chromosome length
		int* numReadsBeg = (int*) calloc(maxChrLen+1, sizeof(int)); // number of reads starting at the position (for preallocation)
		readsGlobal = reads;  // for sorting (comparison)
		
		for(chrI = 0; chrI<nChrom; chrI++)
		{
			// TODO here is some time for a heuristic whether to use counting sort of qsort
			// based on the lenght of the chromosome and the number of reads; for the dense
			// reads we prefer counting sort.
			
			// Counting sort:
			cerr << "Sorting " << chroms[chrI] << endl;
			/** PART 1 of sorting
			 * building an array of lists
			 **/
			int thisChrLen = chrInfo[chroms[chrI]].len;
			int currRead = chrInfo[chroms[chrI]].lastRead;
			int maxNumReads = 0;	// the maximum number of reads at any position
			while(currRead)
			{
				int chosenPosition = (sortMode == '5') ? (reads[currRead].str == '-' ? reads[currRead].end : reads[currRead].beg) : reads[currRead].beg;
				int oldPtr = chromTable[chosenPosition];
				chromTable[chosenPosition] = currRead;
				currRead = reads[currRead].next;		// this and the next line cannot be exchanged
				reads[chromTable[chosenPosition]].next = oldPtr;
				numReadsBeg[chosenPosition] ++;
				if(numReadsBeg[chosenPosition] > maxNumReads)
					maxNumReads ++;	// it can only increase by one.
			}
			/** PART 2 of sorting:
			 * printing and destroying the count list:
			 **/
			if(sortMode == 'b')
			{
				// we want to sort by end position too. for this qsort seems a good solution.
				// allocate the buffer for reads starting at the position
				//cerr << "Allocating " << maxNumReads << endl;
				int* readBuff = (int*) calloc(maxNumReads, sizeof(int));
				for(int pos=0; pos<=thisChrLen; pos++)
				{
					// prepare a table to sort.
					if(chromTable[pos])
					{
						int i=0;
						int foo = numReadsBeg[pos];
						while(chromTable[pos])
						{
							readBuff[i++] = chromTable[pos];
							chromTable[pos] = reads[chromTable[pos]].next;
						}
						//if(numReadsBeg[pos] > 1)	// hack, as usually we will have exactly one read starting there. It is actually slower!!
						//{
							//cerr << "Sorting ends of " << numReadsBeg[pos] << " elements" << endl;
							qsort(readBuff, foo, sizeof(int), cmpReadIndxEnd);
						//}
						// print it out:
						if(fCollapse)
							printf("%s\t%d\t%d\t.\t%g\t+\n", chroms[chrI].c_str(), pos, pos+1, sumList(reads, readBuff, foo));
						else
							for(int j = 0; j < foo; ++j)
								printf("%s", reads[readBuff[j]].line);
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
						if(fCollapse)
							printf("%s\t%d\t%d\t.\t%g\t+\n", chroms[chrI].c_str(), pos, pos+1, sumList2(reads, chromTable, pos));
						else
							while(chromTable[pos])
							{
								printf("%s", reads[chromTable[pos]].line);
								chromTable[pos] = reads[chromTable[pos]].next; 
								// no deallocation here for speed
								//numReadsBeg[pos] = 0;
							}
			}
		}
		time(&end);
		double difftime_sorting = difftime(end, start);
		fprintf(stderr, "Sorting has taken %d seconds\n", (int)(difftime_sorting+0.5));
		return 0;
	}
	catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }
}

int estimateReadNumber(FILE* fh)
{
	// tries to guess how many reads are in the file, so that we
	// can allocate just right the correct amount of RAM. We wont'
	// need reallocating and won't be using too much.
	// This estimate will be biased (underestimated) if the file is 
	// already partly sorted and the top lines are shorter because of
	// lower values
	//return(10);
	fseek(fh, 0, SEEK_END);
	long size = ftell(fh);
	rewind(fh);
	// read the first 1000 lines (or less) and check how many bytes we used
	int nLines = 0;
	int nChars = 0;
	while(nLines<1000)
	{
		int c = fgetc(fh);
		if(c != EOF)
		{
			nChars++;
			if(c == '\n')
			{
				nLines ++;
			}
		}
		else
		{
			break;
		}
	}
	rewind(fh);
	//cout << nLines << "\t" << nChars << "\t" << size << endl;
	int estimate = (int)((double)size * (double)nLines / (double)nChars);
	if(estimate < 10)
	{
		estimate = 10;
	}
	return(estimate);
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

int cmpString(const void* a, const void* b)
{
	return( ((const string*)a)->compare(*(const string*)b) );
}

int cmpSeqreadPtrEnd(const void* a, const void* b)
{
	return( (*(const seqread**)a)->end < (*(const seqread**)b)->end ? -1 : (*(const seqread**)a)->end > (*(const seqread**)b)->end ? 1 : 0 );
}

int cmpReadIndxEnd(const void* a, const void* b)
{
	return ( readsGlobal[*(int*)a].end < readsGlobal[*(int*)b].end ? -1 : readsGlobal[*(int*)a].end > readsGlobal[*(int*)b].end ? 1 : 0) ;
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
