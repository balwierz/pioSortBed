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
| **pioSortBed** | 3.0.0 | `pioSortBed -t 1 input.bed` (single-thread) |
| **pioSortBed** | 3.0.0 | `pioSortBed -t 8 input.bed` (8 threads) |
| **pioSortBed** (low-mem) | 3.0.0 | `pioSortBed --low-mem-ssd -t 1 input.bed` (single-thread) |
| **pioSortBed** (low-mem) | 3.0.0 | `pioSortBed --low-mem-ssd -t 8 input.bed` (8 threads, recommended fast path) |
| **GNU sort** | 9.10 | `LC_ALL=C sort -k1,1 -k2,2n input.bed` (single-thread) |
| **GNU sort** | 9.10 | `LC_ALL=C sort -k1,1 -k2,2n --parallel=8 input.bed` (8 threads) |
| **bedtools** | 2.31.1 | `bedtools sort -i input.bed` |
| **bedops sort-bed** | 2.4.42 | `sort-bed input.bed` |

Wall time and peak RSS (resident set size) measured with GNU time. Times in seconds or milliseconds; memory in MB or GB.


### Wall Time

![Wall time comparison (log-log)](benchmark/benchmark_time.png)

![Wall time comparison (both axes linear)](benchmark/benchmark_time_linear.png)

#### Legend (colour & line style):

Same colour per tool; thread count distinguished by line style (`-t 1` solid, `-t 8` dotted). pioSortBed classic and pioSortBed low-mem are different *algorithms* and get different colours.

| Tool                | Colour     | Marker | Line   | Description |
|---------------------|------------|--------|--------|-------------|
| **pioSortBed 1t**       | <span style="color:#e41a1c">████</span> | ● | solid  | Classic path, single-thread |
| **pioSortBed 8t**       | <span style="color:#e41a1c">████</span> | ● | dotted | Classic path, 8 threads |
| **pioSortBed low-mem 1t** | <span style="color:#c51b7d">████</span> | ◆ | solid  | Low-memory SSD mode, single-thread |
| **pioSortBed low-mem 8t** | <span style="color:#c51b7d">████</span> | ◆ | dotted | Low-memory SSD mode, 8 threads (recommended fast path) |
| **GNU sort 1t**         | <span style="color:#4daf4a">████</span> | ▲ | solid  | Single-thread |
| **GNU sort 8t**         | <span style="color:#4daf4a">████</span> | ▲ | dotted | 8 threads |
| **bedtools**            | <span style="color:#984ea3">████</span> | ✚ | solid  | bedtools sort |
| **bedops**              | <span style="color:#a65628">████</span> | ✦ | solid  | bedops sort-bed |

| Reads | pio 1t ● | pio 8t ● | pio low-mem 1t ◆ | pio low-mem 8t ◆ | sort 1t ▲  | sort 8t ▲  | bedtools ✚ | bedops ✦   |
|------:|---------:|---------:|-----------------:|-----------------:|-----------:|-----------:|-----------:|-----------:|
| 10k   | 0 ms     | 0 ms      | 0 ms             | 0 ms             | 10 ms      | 0 ms       | 10 ms      | 0 ms       |
| 20k   | 0 ms     | 0 ms      | 10 ms            | 10 ms            | 10 ms      | 10 ms      | 20 ms      | 10 ms      |
| 50k   | 10 ms    | 10 ms     | 20 ms            | 10 ms            | 30 ms      | 30 ms      | 40 ms      | 30 ms      |
| 100k  | 10 ms    | 10 ms     | 30 ms            | 20 ms            | 70 ms      | 70 ms      | 90 ms      | 50 ms      |
| 200k  | 30 ms    | 20 ms     | 60 ms            | 30 ms            | 150 ms     | 90 ms      | 170 ms     | 110 ms     |
| 500k  | 100 ms   | 70 ms     | 110 ms           | 70 ms            | 410 ms     | 160 ms     | 450 ms     | 290 ms     |
| 1M    | 210 ms   | 170 ms    | 230 ms           | **150 ms**       | 880 ms     | 320 ms     | 920 ms     | 590 ms     |
| 2M    | 440 ms   | 350 ms    | 550 ms           | **280 ms**       | 1870 ms    | 630 ms     | 1930 ms    | 1360 ms    |
| 5M    | 1120 ms  | 880 ms    | 1670 ms          | **660 ms**       | 5170 ms    | 1670 ms    | 4680 ms    | 3310 ms    |
| 10M   | 2250 ms  | 1780 ms   | 3190 ms          | **1330 ms**      | 11.29 s    | 3450 ms    | 9610 ms    | 6540 ms    |
| 20M   | 4620 ms  | 3640 ms   | 6390 ms          | **2500 ms**      | 24.48 s    | 7120 ms    | 19.93 s    | 14.26 s    |
| 50M   | 19.74 s  | 6150 ms   | 16.34 s          | **5880 ms**      | 1min07.2s  | 19.85 s    | 53.10 s    | 34.42 s    |
| 100M  | 35.66 s  | —         | 32.25 s          | **10.16 s**      | 2min24.0s  | 51.21 s    | —          | 1min07.7s  |
| 200M  | 1min08.4s | —        | 1min09.0s        | **35.16 s**      | 5min28.5s  | 1min45.5s  | —          | 2min34.6s  |

> Sub-50 ms timings (10k–50k) bottom out at GNU `time`'s 10 ms resolution; raw 0 ms readings just mean the tool finished faster than the timer can resolve.
>
> `pio 8t` and `bedtools` are skipped at 100M+ because they would exceed the 30 GB RAM available on this hardware. `pio 8t` allocates a per-chromosome `chromTable` slab on each worker (chr1 alone is ~1 GB), and `bedtools` memory grows linearly with the input. Use `pioSortBed --low-mem-ssd -t 8` instead at those sizes.

**Key observations:**
- **`pioSortBed --low-mem-ssd -t 8` is the fastest configuration at every size from
  1M upwards.** Despite the name, "low-mem" is now also "fast": both passes
  are parallelised, and the per-line index uses a flat 24-byte node table
  instead of a chromosome-length scratch slab. At 200M reads (8.6 GB BED file),
  it's **35.2 s — 3.0× faster than GNU sort 8t, 4.4× over bedops, 9.3× over
  GNU sort 1t**.
- **`pioSortBed --low-mem-ssd -t 1` beats `pioSortBed -t 1` from 50M upwards**
  (50M: 16.3 s vs 19.7 s) and uses ~20% less memory (3.3 GB vs 4.1 GB). At
  smaller sizes the two-pass overhead makes the regular path faster, but the
  low-mem path scales better.
- **`pioSortBed 1t`** (with the LSD radix sort) beats GNU sort 1t by 3–5× across
  the whole range and stays competitive with GNU sort 8t up through 2M.
- **`bedops sort-bed`** remains the closest single-threaded competitor and uses
  the least memory of any tool tested at small sizes.
- The regular `pioSortBed -t 8` parallel bucket-sort path is fast at 50M
  (6.15 s) but its per-thread `chromTable` slabs scale linearly with the
  largest chromosome. For large inputs, `--low-mem-ssd -t 8` strictly dominates.

### Peak Memory (RSS)

![Peak memory usage (log-log)](benchmark/benchmark_memory.png)

![Peak memory usage (both axes linear)](benchmark/benchmark_memory_linear.png)

| Reads | pio 1t ●  | pio 8t ●    | pio low-mem 1t ◆ | pio low-mem 8t ◆ | sort 1t ▲ | sort 8t ▲ | bedtools ✚ | bedops ✦   |
|------:|----------:|------------:|-----------------:|-----------------:|----------:|----------:|-----------:|-----------:|
| 10k   | 6.6 MB    | 6.7 MB      | 6.7 MB           | 7.6 MB           | 3.3 MB    | 3.2 MB    | 8.8 MB     | 2.0 MB     |
| 20k   | 7.7 MB    | 7.1 MB      | 7.3 MB           | 8.5 MB           | 3.2 MB    | 3.3 MB    | 12.4 MB    | 2.5 MB     |
| 50k   | 10.4 MB   | 11.1 MB     | 10.1 MB          | 14.4 MB          | 5.7 MB    | 5.4 MB    | 24.3 MB    | 3.9 MB     |
| 100k  | 15.6 MB   | 16.4 MB     | 15.3 MB          | 18.7 MB          | 9.8 MB    | 9.8 MB    | 44.1 MB    | 6.8 MB     |
| 200k  | 26.1 MB   | 27.0 MB     | 25.9 MB          | 31.8 MB          | 18.5 MB   | 21.3 MB   | 83.8 MB    | 12.0 MB    |
| 500k  | 47.3 MB   | 47.5 MB     | 45.3 MB          | 73.2 MB          | 43.8 MB   | 66.4 MB   | 202.9 MB   | 28.1 MB    |
| 1M    | 90.1 MB   | 82.4 MB     | 77.7 MB          | 121.5 MB         | 86.6 MB   | 161.9 MB  | 400.8 MB   | 54.8 MB    |
| 2M    | 176.3 MB  | 145.8 MB    | 143.4 MB         | 236.3 MB         | 172.7 MB  | 324.1 MB  | 797.8 MB   | 107.9 MB   |
| 5M    | 434.6 MB  | 434.4 MB    | 339.7 MB         | 644.5 MB         | 431.2 MB  | 811.3 MB  | 1.9 GB     | 268.2 MB   |
| 10M   | 865.1 MB  | 865.1 MB    | 667.3 MB         | 1.2 GB           | 861.4 MB  | 1.6 GB    | 3.9 GB     | 535.3 MB   |
| 20M   | 1.7 GB    | 1.7 GB      | 1.3 GB           | 2.5 GB           | 1.7 GB    | 3.2 GB    | 7.8 GB     | 1.0 GB     |
| 50M   | 4.1 GB    | **12.6 GB** | **3.2 GB**       | 5.7 GB           | 4.2 GB    | 8.0 GB    | 19.4 GB    | 2.6 GB     |
| 100M  | 7.2 GB    | —           | **6.5 GB**       | 10.6 GB          | 8.5 GB    | 15.4 GB   | —          | 5.2 GB     |
| 200M  | 13.6 GB   | —           | **13.0 GB**      | 20.5 GB          | 15.4 GB   | 15.4 GB   | —          | 10.4 GB    |

> `pio low-mem 8t` at 100M+ uses `--max-mem=4G` to cap concurrent per-chromosome buffers.

**Key observations:**
- **`pioSortBed --low-mem-ssd -t 8` is the recommended fast path for files.**
  At 50M+ it beats `pio -t 8` on both wall time and memory (50M: 5.88 s vs
  6.15 s, 5.7 GB vs 12.6 GB). At 100M and 200M it's the only pioSortBed mode
  besides single-threaded that fits in 32 GB RAM at all.
- **`pioSortBed --low-mem-ssd -t 1` is the lowest-memory pioSortBed mode** —
  matches or beats `pio -t 1` on memory at every size from 5M up (50M: 3.2 GB
  vs 4.1 GB) and at 200M reads is the lowest-RAM pioSortBed mode at 13.0 GB.
- **Memory grows ~linearly with input for every tool.** At 200M: pio-lm 8t
  21 GB, pio-lm 1t 13 GB, GNU sort 15 GB, bedops 10 GB. `bedtools` and
  `pio -t 8` would have needed >32 GB and were skipped.
- **`bedops sort-bed`** uses the least memory throughout — a sensible choice
  on RAM-constrained systems where wall time isn't critical.
- **`pioSortBed -t 8`** trades RAM for speed via per-thread `chromTable`
  slabs. The `--max-mem=N[GMK]` flag (v2.2.0+) lets you cap the peak in
  exchange for wall time. For large inputs `--low-mem-ssd -t 8` is now
  strictly better and `pio -t 8` is mainly useful below ~50M.

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

> The `pioSortBed low-mem` row was measured with the default thread count (= all cores on the bench box) while the other "8t" rows used `-t 8` / `--parallel=8`. The synthetic table above splits low-mem into `-t 1` and `-t 8`; this real-data table will be re-run on the next benchmark cycle.

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

> Same caveat as the chr20 table above: the `pioSortBed low-mem` row used the default thread count, not `-t 8`. Both real-data tables will be re-run with explicit `-t 1` / `-t 8` low-mem rows on the next benchmark cycle.

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
