# pioSortBed

**Ultra-fast BED file sorter for genomics**

Sorts BED files by chromosome and start coordinate, equivalent to:
```
LC_ALL=C sort -k1,1 -k2,2n file.bed
```
but significantly faster on large datasets. Supports BED3, BED6, and extended BED formats.

## Algorithm

pioSortBed uses a hybrid sort strategy:
- **Files with < 10M reads** (configurable via `--bucket-cutoff`): classic **O(n log n)** comparison sort (`std::sort` on an index array). Optionally parallel via `--threads`.
- **Files with ≥ 10M reads**: bucket sort (counting sort), which avoids coordinate comparisons entirely — reads are placed directly into position-indexed buckets. This gives **O(n + m)** complexity, where *n* is the number of regions and *m* is the maximum chromosome length.

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
| `--bucket-cutoff N` | Use bucket sort for files with ≥N reads (default: 10M; 0 = always bucket sort) |
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

Sorting random BED6 files (10 chromosomes, coordinates 0–249 Mbp) on a laptop with 32 GB RAM. Wall time and peak RSS measured with GNU time. All tools verified to produce identical sort order.

![Benchmark plot](benchmark_plot.png)

Solid lines = wall time (left axis), dashed lines = peak memory (right axis). Log–log scale.

### Wall time

| Reads | pioSortBed v2 | GNU sort | bedtools sort | bedops sort-bed | pioSortBed v1 |
|------:|--------------:|---------:|--------------:|----------------:|--------------:|
| 10k   | **< 10 ms**   | 10 ms    | 10 ms         | 10 ms           | 3.96 s        |
| 100k  | **30 ms**     | 130 ms   | 180 ms        | 110 ms          | 3.70 s        |
| 1M    | **0.41 s**    | 0.53 s   | 1.59 s        | 1.10 s          | 5.28 s        |
| 5M    | **2.44 s**    | 2.80 s   | 8.26 s        | 5.78 s          | 11.11 s       |
| 10M   | 12.31 s       | **5.83 s** | 16.38 s     | 11.47 s         | 19.66 s       |

### Peak memory (RSS)

| Reads | pioSortBed v2 | GNU sort | bedtools sort | bedops sort-bed | pioSortBed v1 |
|------:|--------------:|---------:|--------------:|----------------:|--------------:|
| 10k   | **3.3 MB**    | 4.3 MB   | 8.6 MB        | 1.9 MB          | 1.9 GB        |
| 100k  | **9.9 MB**    | 12.0 MB  | 48.3 MB       | 8.6 MB          | 1.9 GB        |
| 1M    | **72.9 MB**   | 164.7 MB | 415.0 MB      | 56.7 MB         | 1.9 GB        |
| 5M    | **356.2 MB**  | 814.3 MB | 2.0 GB        | 268.1 MB        | 2.3 GB        |
| 10M   | 1.5 GB        | 1.6 GB   | 3.9 GB        | **535.4 MB**    | 2.7 GB        |

> **pioSortBed v1** always allocates a ~2 GB bucket table regardless of input size. **v2** uses a hybrid strategy (classic sort for < 10M reads, bucket sort above), giving both better speed and lower memory for typical files.

Run `bash benchmark.sh` to reproduce (requires GNU time, gnuplot for the plot).

## Compile-time Limits

These constants can be changed and the program recompiled if needed:

| Constant | Default | Description |
|----------|---------|-------------|
| `lineBufSize` | 1024 bytes | Maximum BED line length (stdin only; no limit for file input) |
| `chrNameBufSize` | 256 bytes | Maximum chromosome name length |
| `chrLenLimit` | 1 Gbp | Maximum chromosome/contig length |
| `defaultBucketCutoff` | 10M | Hybrid sort threshold (overridden by `--bucket-cutoff`) |

## Author

Piotr Balwierz — Imperial College London
