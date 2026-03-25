# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make
```

No external dependencies — CLI11 (header-only argument parser) is included in the repository. The Makefile links statically by default. To build without static linking, remove `-static` from `LDFLAGS` in the Makefile.

Manual compilation:
```bash
g++ pioSortBed.cpp -o pioSortBed -O3 -std=c++11 -fopenmp -static
```

## Usage

```bash
./pioSortBed file.bed > sorted.bed
./pioSortBed - < file.bed > sorted.bed   # read from stdin
```

Options: `--sort s` (by start, default), `--sort b` (by start+end), `--sort 5` (by 5' end respecting strand), `--collapse` (sum weights of overlapping regions), `--ral` (RAL input format instead of BED), `--bucket-cutoff N` (use bucket sort for files with ≥N reads; default 50M; 0 = always bucket sort).

## Architecture

Single-file C++ program (`pioSortBed.cpp`). Hybrid sort strategy: files with fewer than 50M reads (configurable via `--bucket-cutoff`) use classic O(n log n) comparison sort (`std::sort` on an index array), while larger files use bucket/counting sort giving O(n+m) complexity where m = max chromosome length.

**Data flow:**
1. Parse input BED/RAL lines into `seqread` structs (stored in a flat array, dynamically grown with `realloc`)
2. Per chromosome, maintain a linked list of read indices (`seqread.next` field, head stored in `chrInfoT.lastRead`)
3. Sort chromosomes alphabetically (`std::sort` on chromosome name strings)
4. **Classic sort path** (< cutoff reads): Walk linked lists to stamp integer chromosome indices (`chrIdx`, sharing storage with `next` via a union), collect all read indices into a flat array, `std::sort` with inlined lambda comparator, linear print.
5. **Bucket sort path** (≥ cutoff reads): For each chromosome, scatter reads into `chromTable` array (indexed by position), walk from 0 to chromosome length printing reads in order.

**Key compile-time limits** (defined as constants at top of file):
- `lineBufSize` = 1024 — max BED line length (stdin only; no limit for mmap)
- `chrNameBufSize` = 256 — max chromosome name length
- `chrLenLimit` = 1,000,000,000 — max coordinate value (1 Gbp); increase and recompile for non-standard genomes
- `defaultBucketCutoff` = 50,000,000 — hybrid sort threshold (overridden by `--bucket-cutoff`)

**Memory:** Stores all data in RAM. Classic sort path uses ~readCount×4 bytes extra. Bucket sort path allocates chromTable[maxChrLen+1] (up to 4 GB).
