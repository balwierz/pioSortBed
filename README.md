# pioSortBed

**Ultra-fast BED file sorter for genomics**

Sorts BED files by chromosome and start coordinate, equivalent to:
```
LC_ALL=C sort -k1,1 -k2,2n file.bed
```
but significantly faster on large datasets. Supports BED3, BED6, and extended BED formats.

## Algorithm

pioSortBed uses a hybrid sort strategy:
- **Files with < 50M reads** (configurable via `--bucket-cutoff`): classic **O(n log n)** comparison sort (`std::sort` on an index array). Optionally parallel via `--threads`.
- **Files with ≥ 50M reads**: bucket sort (counting sort), which avoids coordinate comparisons entirely — reads are placed directly into position-indexed buckets. This gives **O(n + m)** complexity, where *n* is the number of regions and *m* is the maximum chromosome length.

> Bucket sort has overhead proportional to chromosome length (up to 4 GB allocation), so the classic sort path is preferred for smaller files.

## Installation

**Dependencies:** GCC (no external libraries needed — CLI11 is bundled)

```bash
make
```

Manual compilation:
```bash
g++ pioSortBed.cpp -o pioSortBed -O3 -std=c++11 -fopenmp -static
```

## Usage

```
pioSortBed [options] <input.bed>
pioSortBed [options] -   # read from standard input
```

| Option | Description |
|--------|-------------|
| `-s s` / `--sort s` | Sort by start coordinate (default) |
| `-s b` / `--sort b` | Sort by start and end coordinate |
| `-s 5` / `--sort 5` | Sort by 5' end (respects strand: col 6) |
| `-r` / `--ral` | Input is in RAL format instead of BED |
| `--collapse` | Collapse overlapping regions, summing weights |
| `--bucket-cutoff N` | Use bucket sort for files with ≥N reads (default: 50M; 0 = always bucket sort) |
| `-t N` / `--threads N` | Number of threads for classic sort (0 = all cores; 1 = single-threaded) |
| `-h` / `--help` | Show help message |

**Examples:**
```bash
pioSortBed input.bed > sorted.bed
cat input.bed | pioSortBed - > sorted.bed
pioSortBed --sort b input.bed > sorted.bed
```

## Memory Requirements

All data is loaded into memory. Expect approximately **2× the input file size** in RAM usage.

## Benchmark

Sorting random BED6 files (10 chromosomes, coordinates 0–249 Mbp). Wall time and peak RSS measured with GNU time. All tools verified to produce identical sort order.

### Wall time

![Wall time](benchmark_time.png)

| Reads | pio 1t | pio 8t | sort 1t | sort 8t | bedtools | bedops |
|------:|-------:|-------:|--------:|--------:|---------:|-------:|
| 10k   | < 10 ms | < 10 ms | 30 ms | 30 ms | 10 ms | 10 ms |
| 100k  | **30 ms** | **30 ms** | 170 ms | 170 ms | 180 ms | 110 ms |
| 1M    | 490 ms | **420 ms** | 2.02 s | 1.17 s | 1.79 s | 1.23 s |
| 5M    | 2.91 s | **2.36 s** | 12.15 s | 7.20 s | 9.45 s | 6.33 s |
| 10M   | 6.02 s | **5.02 s** | 26.10 s | 15.37 s | 19.62 s | 12.92 s |
| 50M   | **26.92 s** | 27.02 s | 2min33s | 1min27s | 1min47s | 1min06s |
| 100M  | **51.81 s** | 52.25 s | 5min36s | 3min14s | 3min46s | 2min16s |

### Peak memory (RSS)

![Peak memory](benchmark_memory.png)

| Reads | pio 1t | pio 8t | sort 1t | sort 8t | bedtools | bedops |
|------:|-------:|-------:|--------:|--------:|---------:|-------:|
| 10k   | 3.5 MB | 3.7 MB | 5.4 MB | 3.6 MB | 7.2 MB | **2.7 MB** |
| 100k  | **10.0 MB** | 10.1 MB | 14.2 MB | 12.9 MB | 47.1 MB | 11.2 MB |
| 1M    | 69.9 MB | 72.2 MB | 89.2 MB | 165.7 MB | 410.1 MB | **55.5 MB** |
| 5M    | 339.1 MB | 355.4 MB | 433.7 MB | 818.2 MB | 1.9 GB | **269.2 MB** |
| 10M   | 674.4 MB | 711.1 MB | 864.5 MB | 1.6 GB | 3.9 GB | **536.2 MB** |
| 50M   | 4.1 GB | 4.1 GB | 4.2 GB | 8.0 GB | 19.4 GB | **2.6 GB** |
| 100M  | 7.2 GB | 7.2 GB | 8.5 GB | 15.9 GB | 38.7 GB | **5.2 GB** |

> pioSortBed uses a hybrid strategy: classic sort for < 50M reads, bucket sort above (configurable via `--bucket-cutoff`). At 50M reads the bucket sort kicks in — note the flat time scaling and constant memory per read. At 100M reads pioSortBed is 4–6× faster than alternatives. GNU sort 8-thread uses 2× the memory of its single-threaded mode.

Run `bash benchmark.sh` to reproduce (requires GNU time; gnuplot for the plot).

## Compile-time Limits

These constants can be changed and the program recompiled if needed:

| Constant | Default | Description |
|----------|---------|-------------|
| `lineBufSize` | 1024 bytes | Maximum BED line length (stdin only; no limit for file input) |
| `chrNameBufSize` | 256 bytes | Maximum chromosome name length |
| `chrLenLimit` | 1 Gbp | Maximum chromosome/contig length |
| `defaultBucketCutoff` | 50M | Hybrid sort threshold (overridden by `--bucket-cutoff`) |

## Author

Piotr Balwierz
