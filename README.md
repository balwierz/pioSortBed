# pioSortBed

**Ultra-fast BED file sorter for genomics**

Sorts BED files by chromosome and start coordinate, equivalent to:
```
LC_ALL=C sort -k1,1 -k2,2n file.bed
```
but significantly faster on large datasets. Supports BED3, BED6, and extended BED formats.

## Algorithm

pioSortBed uses a hybrid sort strategy:
- **Files with < 50M reads** (configurable via `--bucket-cutoff`): classic **O(n log n)** comparison sort (`std::sort` on an index array). Single-threaded by default; parallel via TBB-backed `std::sort(std::execution::par, ...)` at `--threads > 1`.
- **Files with ≥ 50M reads**: bucket sort (counting sort), which avoids coordinate comparisons entirely — reads are placed directly into position-indexed buckets. This gives **O(n + m)** complexity, where *n* is the number of regions and *m* is the maximum chromosome length.

Both phases parallelize:
- **Parsing** (mmap path): file is split into newline-aligned chunks parsed concurrently into per-chunk per-chromosome partial linked lists, merged in a final serial pass. Header lines are emitted to stdout in a leading pre-pass before chunking.
- **Bucket sort** (multi-thread): chromosomes processed concurrently, each thread allocating its own per-chromosome `chromTable` sized exactly to that chromosome's max coordinate. Output goes through a producer-consumer barrier so per-chromosome buffers are flushed to stdout in alphabetical order as soon as their turn comes up — capping live output buffers to roughly the worker count.

> Bucket sort has overhead proportional to chromosome length, and the parallel bucket-sort path also pays for per-thread `chromTable` slabs and queued output buffers. Use `--threads 1` if RAM is tight; the single-threaded path keeps a flat ~4 GB ceiling at 50M reads.

## Installation

**Dependencies:** GCC ≥ 9 (C++17), oneTBB (`libtbb-dev` on Debian/Ubuntu).
CLI11 is bundled in this repo.

```bash
make
make test       # optional: runs the test suite
make install    # optional: installs to /usr/local/bin (override PREFIX=...)
```

Manual compilation:
```bash
g++ src/pioSortBed.cpp -Isrc -o pioSortBed -O3 -std=c++17 \
    -static-libstdc++ -static-libgcc -ltbb -DVERSION_STRING=\"2.1.0\"
```

> Parallelism uses C++17 `std::execution::par` algorithms backed by oneTBB,
> not OpenMP. The C++ runtime is linked statically; libtbb stays dynamic.

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
| `-n` / `--natural-sort` | Natural chromosome order: `chr2 < chr10` (default: lexicographic) |
| `-r` / `--ral` | Input is in RAL format instead of BED |
| `--collapse` | Collapse overlapping regions, summing weights |
| `--low-mem-ssd` | Low-memory two-pass file mode (SSD-friendly). Slower than default, but lower peak RAM. Requires file input (not stdin or gzip). |
| `--bucket-cutoff N` | Use bucket sort for files with ≥N reads (default: 50M; 0 = always bucket sort) |
| `-t N` / `--threads N` | Number of threads for classic sort (0 = all cores; 1 = single-threaded) |
| `--max-mem=N[GMK]` | Memory budget for the parallel bucket-sort path (e.g. `4G`, `500M`). Caps concurrent per-chromosome `chromTable` allocations so peak RAM stays within budget. Default: no cap. |
| `-h` / `--help` | Show help message |

BED header lines (`track`, `browser`, `#` comments) are passed through unchanged to output. Gzip-compressed input (`.gz` extension) is transparently decompressed.

**Examples:**
```bash
pioSortBed input.bed > sorted.bed
pioSortBed input.bed.gz > sorted.bed          # gzip input
cat input.bed | pioSortBed - > sorted.bed
pioSortBed --sort b input.bed > sorted.bed
pioSortBed --natural-sort input.bed > sorted.bed   # chr2 before chr10
```

## Memory Requirements

All data is loaded into memory. Expect approximately **2× the input file size** in RAM usage.

## Benchmark Results

Comprehensive sorting benchmark on realistic BED6 files (10 chromosomes, coordinates 0–249 Mbp). All tools verified to produce identical sort order.

### System Configuration

**Hardware:** Lenovo ThinkPad P1 Gen7
- CPU: Intel Core Ultra 7 155H (Meteor Lake, Intel 4) — 16 cores / 22 threads: 6 P-cores (Redwood Cove, up to 4.8 GHz) + 8 E-cores (Crestmont, up to 3.8 GHz) + 2 LP E-cores (Crestmont, up to 2.5 GHz)
- RAM: 32 GB LPCAMM2 @ 7467 MT/s
- SSD: KIOXIA KXG8AZNV1T02 NVMe — random 4 kB read: 177 MiB/s, 45.3k IOPS (fio: `--rw=randread --bs=4k --size=1G --numjobs=4 --runtime=30`)

**Tools & Command Lines:**

| Tool | Version | Command |
|------|---------|---------|
| **pioSortBed** | 2.2.1 | `pioSortBed -t 1 input.bed` (single-thread) |
| **pioSortBed** | 2.2.1 | `pioSortBed -t 8 input.bed` (8 threads) |
| **pioSortBed** (low-mem) | 2.2.1 | `pioSortBed --low-mem-ssd input.bed` (two-pass SSD-friendly mode) |
| **GNU sort** | 9.10 | `LC_ALL=C sort -k1,1 -k2,2n input.bed` (single-thread) |
| **GNU sort** | 9.10 | `LC_ALL=C sort -k1,1 -k2,2n --parallel=8 input.bed` (8 threads) |
| **bedtools** | 2.31.1 | `bedtools sort -i input.bed` |
| **bedops sort-bed** | 2.4.42 | `sort-bed input.bed` |

Wall time and peak RSS (resident set size) measured with GNU time. Times in seconds or milliseconds; memory in MB or GB.


### Wall Time

![Wall time comparison](benchmark/benchmark_time.png)

![Wall time comparison (linear scale)](benchmark/benchmark_time_linear.png)

#### Legend (point style & color):

| Tool                | Color      | Point Style | Description |
|---------------------|------------|-------------|-------------|
| **pioSortBed 1t**   | <span style="color:#e41a1c">████</span> | ●           | Single-thread (default) |
| **pioSortBed 8t**   | <span style="color:#377eb8">████</span> | ■           | 8 threads |
| **pioSortBed low-mem** | <span style="color:#f781bf">████</span> | ◆        | Low-memory SSD mode |
| **GNU sort 1t**     | <span style="color:#4daf4a">████</span> | ▲           | Single-thread |
| **GNU sort 8t**     | <span style="color:#ff7f00">████</span> | ▼           | 8 threads |
| **bedtools**        | <span style="color:#984ea3">████</span> | ✚           | bedtools sort |
| **bedops**          | <span style="color:#a65628">████</span> | ✦           | bedops sort-bed |

| Reads | pio 1t ● | pio 8t ■   | pio low-mem ◆ | sort 1t ▲  | sort 8t ▼ | bedtools ✚ | bedops ✦ |
|------:|----------|------------|---------------|------------|-----------|------------|-----------|
| 10k   | 0 ms     | 0 ms       | 0 ms          | 0 ms       | 10 ms     | 10 ms      | 0 ms      |
| 20k   | 0 ms     | 0 ms       | 0 ms          | 10 ms      | 10 ms     | 10 ms      | 10 ms     |
| 50k   | 10 ms    | 0 ms       | 10 ms         | 40 ms      | 30 ms     | 40 ms      | 30 ms     |
| 100k  | 20 ms    | 10 ms      | 10 ms         | 90 ms      | 70 ms     | 90 ms      | 60 ms     |
| 200k  | 30 ms    | 40 ms      | 30 ms         | 170 ms     | 100 ms    | 170 ms     | 120 ms    |
| 500k  | 130 ms   | 70 ms      | 100 ms        | 410 ms     | 170 ms    | 440 ms     | 290 ms    |
| 1M    | 210 ms   | 170 ms     | 200 ms        | 860 ms     | 290 ms    | 910 ms     | 600 ms    |
| 2M    | 430 ms   | 360 ms     | 400 ms        | 1880 ms    | 620 ms    | 1880 ms    | 1280 ms   |
| 5M    | 1150 ms  | 910 ms     | 1010 ms       | 5200 ms    | 1650 ms   | 4680 ms    | 3290 ms   |
| 10M   | 2280 ms  | 1780 ms    | 2010 ms       | 11.25 s    | 3440 ms   | 9580 ms    | 6590 ms   |
| 20M   | 4600 ms  | 3620 ms    | 4240 ms       | 24.43 s    | 7290 ms   | 19.78 s    | 14.10 s   |
| 50M   | 18.95 s  | **6.15 s** | 10.62 s       | 1min07.4s  | 19.78 s   | 54.25 s    | 34.33 s   |

> Sub-50ms timings (10k–50k) bottom out at GNU `time`'s 10 ms resolution; raw 0 ms readings just mean the tool finished faster than the timer can resolve.

**Key observations:**
- **pioSortBed 8t** is fastest from 50k upwards. The parallel bucket-sort
  at 50M shines: **3.1× over pio 1t, 3.2× over GNU sort 8t, 5.5× over
  bedops, 8.5× over bedtools, 10.8× over GNU sort 1t**.
- **pioSortBed 1t** (with the v2.1.1 LSD radix sort + v2.2.1 perf polish) now
  beats GNU sort 1t by **3–5×** across the entire range and stays
  competitive with GNU sort 8t through 2M.
- **pioSortBed low-mem mode** is comparable to `-t 8` at 10M–20M and uses
  the least RAM of the pio variants — useful when `-t 8`'s memory peak
  matters more than the last 10–15% of speed.
- **bedops sort-bed** remains the closest single-threaded competitor and
  uses the lowest memory of any tool tested.
- **GNU sort 8t** scales well at 1M–10M but cannot match per-chromosome
  parallelism on a single-file sort at 50M.

### Peak Memory (RSS)

![Peak memory usage](benchmark/benchmark_memory.png)

![Peak memory usage (linear scale)](benchmark/benchmark_memory_linear.png)

| Reads | pio 1t ●  | pio 8t ■    | pio low-mem ◆ | sort 1t ▲ | sort 8t ▼ | bedtools ✚ | bedops ✦ |
|------:|-----------|-------------|---------------|-----------|-----------|------------|-----------|
| 10k   | 6.7 MB    | 6.7 MB      | 6.9 MB        | 3.3 MB    | 3.3 MB    | 8.9 MB     | 2.0 MB    |
| 20k   | 7.8 MB    | 8.1 MB      | 7.6 MB        | 3.4 MB    | 3.1 MB    | 12.6 MB    | 2.5 MB    |
| 50k   | 10.7 MB   | 11.2 MB     | 10.6 MB       | 5.7 MB    | 5.8 MB    | 24.5 MB    | 4.1 MB    |
| 100k  | 15.6 MB   | 16.2 MB     | 15.4 MB       | 9.8 MB    | 9.8 MB    | 44.2 MB    | 6.6 MB    |
| 200k  | 26.0 MB   | 27.1 MB     | 25.9 MB       | 18.1 MB   | 21.1 MB   | 83.8 MB    | 12.1 MB   |
| 500k  | 47.5 MB   | 47.8 MB     | 44.8 MB       | 43.9 MB   | 66.3 MB   | 203.1 MB   | 28.1 MB   |
| 1M    | 90.1 MB   | 83.0 MB     | 77.3 MB       | 86.7 MB   | 161.8 MB  | 401.0 MB   | 54.8 MB   |
| 2M    | 176.4 MB  | 146.1 MB    | 138.8 MB      | 172.6 MB  | 323.9 MB  | 798.0 MB   | 108.0 MB  |
| 5M    | 434.6 MB  | 433.9 MB    | 316.9 MB      | 430.8 MB  | 810.9 MB  | 1.9 GB     | 268.1 MB  |
| 10M   | 865.0 MB  | 865.0 MB    | 621.5 MB      | 861.4 MB  | 1.6 GB    | 3.9 GB     | 535.3 MB  |
| 20M   | 1.7 GB    | 1.7 GB      | 1.2 GB        | 1.7 GB    | 3.2 GB    | 7.8 GB     | 1.0 GB    |
| 50M   | 4.1 GB    | **12.6 GB** | 3.0 GB        | 4.2 GB    | 8.0 GB    | 19.4 GB    | 2.6 GB    |

**Key observations:**
- **pioSortBed 8t at 50M (12.6 GB)** is the price for the 3.1× speedup —
  each thread allocates its own per-chromosome `chromTable` and queued
  per-chromosome output buffers accumulate. Below 50M, 8t and 1t memory
  are within a few %; the classic-sort path doesn't allocate per-chromosome
  scratch. The new `--max-mem=N[GMK]` flag (v2.2.0+) lets you cap this
  peak at the cost of some wall time — see CHANGELOG for the trade-off table.
- **pioSortBed low-mem mode** keeps a flat ~3 GB ceiling at 50M by
  processing one chromosome at a time — pick this when both wall time and
  RAM matter and you're on SSD.
- **bedops** achieves the lowest memory on small-to-medium files and is
  competitive throughout.
- **GNU sort 8t** uses ~2× the RAM of single-threaded sort (thread-local buffers).
- **bedtools** memory grows ~linearly with input — 19.4 GB at 50M.

### Performance Summary

**pioSortBed strategy:**
- **Files < 50M reads** (configurable `--bucket-cutoff`): Classic **O(n log n)** comparison sort on an index array
  - Single-threaded `std::sort` with inlined comparator at `-t 1`
  - Parallel `std::sort(std::execution::par, ...)` (TBB) at `-t > 1`
- **Files ≥ 50M reads**: Bucket/counting sort — O(n + m) complexity where m = max chromosome length
  - Single-thread: one shared `chromTable[maxChrLen+1]` reused across chromosomes
  - Multi-thread (v2.1.0): chromosomes processed concurrently, each thread allocates its own `chromTable[thisChrLen+1]`. Per-chrom output is buffered through `open_memstream` and flushed under a producer-consumer barrier so output is written in alphabetical order
- **Parallel parsing** (mmap path): file split into newline-aligned chunks parsed in parallel; per-chunk per-chromosome partials merged at the end
- **Low-memory mode** (`--low-mem-ssd`): Two-pass algorithm for RAM-constrained environments
  - Pass 1: Scan file once, store line offsets per chromosome (minimal RAM)
  - Pass 2: Process one chromosome at a time, sorting and printing in isolation
  - Trade-off: slower than default mode, but peak RAM ∝ largest chromosome (not whole file)
  - Best for large genomic files on small-RAM systems

To reproduce: `bash benchmark/benchmark.sh` (requires GNU time; gnuplot for plots).

### Real-data Benchmark: NA12878 WGS (chr20, 120M reads)

Benchmark on real Illumina WGS reads: NA12878 (HG001) 300x HiSeq, chr20, aligned to GRCh38 (GIAB/NHGRI). 120,499,538 reads, 7.9 GB BED file. All tools use the bucket-sort path (>50M reads).

| Tool | Wall time | Peak RSS |
|------|-----------|----------|
| **pioSortBed 1t** | 9.4 s | 10.8 GB |
| **pioSortBed 8t** | 9.2 s | 10.8 GB |
| **pioSortBed low-mem** | 12.0 s | 15.5 GB |
| **GNU sort 1t** | 1min 09.8s | 13.2 GB |
| **GNU sort 8t** | 33.9 s | 22.2 GB |
| **bedops sort-bed** | 59.3 s | 9.9 GB |
| **bedtools sort** | 6min 14.3s | 40.3 GB |

**pioSortBed is 7.4× faster than GNU sort (single-thread) and 3.6× faster than GNU sort (8-thread).** bedops is competitive on memory but 6.3× slower than pioSortBed. bedtools is the slowest and uses the most RAM (40.3 GB).

To reproduce: `bash benchmark/benchmark_na12878.sh` (streams ~12 GB from NCBI FTP on first run).

### Real-data Benchmark: NA12878 WGS (all chromosomes, exactly 100M reads)

100,000,000 reads randomly sampled from all standard chromosomes (chr1–22, X, Y) of the same HG001 GRCh38 300x BAM. Sampled by streaming a 2% subsample (~114M reads) then `shuf -n 100000000`. Reads span all chromosomes — realistic multi-chromosome sort workload.

| Tool | Wall time | Peak RSS |
|------|-----------|----------|
| **pioSortBed 1t** | 38.0 s | 9.7 GB |
| **pioSortBed 8t** | 37.2 s | 9.7 GB |
| **pioSortBed low-mem** | 39.7 s | 8.6 GB |
| **GNU sort 1t** | 3min 07.5s | 11.0 GB |
| **GNU sort 8t** | 57.0 s | 17.3 GB |
| **bedops sort-bed** | 1min 10.2s | 8.2 GB |
| **bedtools sort** | 6min 24.7s | 40.1 GB |

**pioSortBed is 4.9× faster than GNU sort (single-thread) and 1.5× faster than GNU sort (8-thread).** bedops has the lowest memory (8.2 GB) but is 1.9× slower. bedtools is slowest and uses the most RAM.

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
