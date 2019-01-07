Fastest implementation of BED file sorting [genomics]

It replaces the UNIX sort -k1,1 -k2,2n file.bed (in LC_ALL=C)
and the bedtools sort commands.

It uses (what I believe is called) bucket sorting. So effectively does not do any coordinate comparisons,
just the indexing. This way it needs to store all data in the memory. Roughly twice as much memory needed than
the size of the BED file. There are some compile-time limits on the lenghts of lines [1024] and chromosome
name lengths [256], and chromosome length limits [1Gbp]. You can change these and recompile.

It can do some trivial operations like collapsing regions if they are multiple lines regions with the same coordinates.

It probably is compatible with Unicode characters in read names :-)


Piotr Balwierz
Imperial College London
