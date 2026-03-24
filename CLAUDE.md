# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make
```

Requires Boost libraries (`boost_program_options`). The Makefile links statically by default. To build without static linking, remove `-static` from `LDFLAGS` in the Makefile.

If Boost is in a non-standard location, pass include/lib paths:
```bash
g++ -I/path/to/boost/include -L/path/to/boost/lib pioSortBed.cpp -o pioSortBed -O3 -lboost_program_options -static
```

## Usage

```bash
./pioSortBed file.bed > sorted.bed
./pioSortBed - < file.bed > sorted.bed   # read from stdin
```

Options: `--sort s` (by start, default), `--sort b` (by start+end), `--sort 5` (by 5' end respecting strand), `--collapse` (sum weights of overlapping regions), `--ral` (RAL input format instead of BED).

## Architecture

Single-file C++ program (`pioSortBed.cpp`). The algorithm is bucket sort (counting sort), giving O(n+m) complexity where n = number of reads and m = max chromosome length — faster than `sort -k1,1 -k2,2n` for large datasets.

**Data flow:**
1. Parse input BED/RAL lines into `seqread` structs (stored in a flat array, dynamically grown with `realloc`)
2. Per chromosome, maintain a linked list of read indices (`seqread.next` field, head stored in `chrInfoT.lastRead`)
3. Sort chromosomes alphabetically (`std::sort` on chromosome name strings)
4. For each chromosome: scatter reads into `chromTable` array (indexed by position) — this is the bucket sort step
5. Walk `chromTable` from 0 to chromosome length, printing reads in order

**Key compile-time limits** (defined as constants at top of file):
- `lineBufSize` = 1024 — max BED line length
- `chrNameBufSize` = 256 — max chromosome name length
- `chrLenLimit` = 1,000,000,000 — max coordinate value (1 Gbp); increase and recompile for non-standard genomes

**Memory:** Stores all data in RAM. Roughly 2× the BED file size is required.

**Global:** `readsGlobal` is a global pointer used as a workaround to pass data into the `qsort` comparator `cmpReadIndxEnd` (used only in `--sort b` mode).
