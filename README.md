Fastest implementation of BED file sorting [genomics]

It replaces the UNIX sort -k1,1 -k2,2n file.bed (in LC_ALL=C)
and the bedtools sort commands.

Input files: BED3, BED6+n etc

It uses (what I believe is called) bucket sorting. So effectively does not do any coordinate comparisons,
just the indexing. In this way it is not anymore O(n*log(n)) problem, but O(n+m) problem, where n is the numer of reads
and m it the maximum length of a chromosome/contig. So for large datasets it runs in linear time of the number of reads.
But it sucks on really small files


pioSortBed needs to *store all data in the memory*. Roughly twice as much memory needed than
the size of the BED file. It is possible to add some swapping in the future.

There are some compilation-time limits on the lenghts of lines [1024], chromosome name lengths [256]
and chromosome length limits [1Gbp]. You can change these and recompile.

It can do some trivial operations like collapsing regions if they are multiple lines regions with the same coordinates.
It does not aim at replacing bedops, bedtools, GenomicRanges etc

It is probably compatible with Unicode characters in read names :-) Uses Boost for hast tables and command line opions.


Piotr Balwierz
Imperial College London
