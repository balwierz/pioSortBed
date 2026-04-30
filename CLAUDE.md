# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make            # produces ./pioSortBed (links libtbb.so dynamically)
make test       # runs the 29-check test suite under test/test.sh
make install    # installs to $PREFIX/bin (default /usr/local)
```

**Dependencies:** GCC ≥ 9 (C++17), oneTBB (`libtbb-dev` on Debian/Ubuntu).
CLI11 is bundled in `src/CLI11.hpp`.

The default build links libstdc++ and libgcc statically (`-static-libstdc++ -static-libgcc`) but leaves libtbb dynamic — Debian/Ubuntu only ships `libtbb.so`, not `libtbb.a`. For the GitHub release artefact (a fully self-contained binary) build a static oneTBB locally and use `make release-binary TBB_LIB=/path/to/libtbb.a`. The Makefile comments include the oneTBB cmake recipe.

Manual compilation:
```bash
g++ src/pioSortBed.cpp -Isrc -o pioSortBed -O3 -std=c++17 \
    -static-libstdc++ -static-libgcc -ltbb -DVERSION_STRING=\"3.0.12\"
```

## Usage

```bash
./pioSortBed file.bed > sorted.bed
./pioSortBed file.bed.gz > sorted.bed     # gzip transparently decompressed (slurped to memory)
./pioSortBed - < file.bed > sorted.bed    # stdin (also slurped)
./pioSortBed -o sorted.bed file.bed       # explicit output file
```

Options:
- `--sort s` (by start, default), `--sort b` (start + end), `--sort 5` (5'-end, strand-aware)
- `-r` / `--ral` — RAL input format instead of BED
- `-n` / `--natural-sort` — `chr2 < chr10` instead of lexicographic
- `--collapse` — sum weights of regions sharing (chr, start)
- `--low-mem-ssd` — two-pass mode, **recommended fast path** for any input ≥ ~1 M reads
- `--bucket-cutoff N` — opt in to bucket sort within the classic path (`-1` = disabled, default; `0` = always; `N ≥ 1` = at ≥ N reads)
- `-t N` / `--threads N` — worker pool size (`0` = all cores, default; `1` = serial)
- `--max-mem=N[GMK]` — per-chromosome scratch cap (safety, not optimisation; see README's `--max-mem` sweep section)
- `-o FILE` / `--output FILE` — write to file (default: stdout)
- `-v` / `--verbose` — opt in to timing / chromosome-length info on stderr (default: silent)

BED header lines (`#`, `track `, `browser `) pass through unchanged.

## Architecture

Single source file (`src/pioSortBed.cpp`, ~2700 lines). Section banners divide it into:
**CORE TYPES & ALLOCATOR** / **PARSERS** / **LOW-MEMORY SSD PATH** / **CLASSIC SORT PATH** / **Parallel mmap parsing** / **CLI / MAIN**. Grep `^// =\+$` to jump between sections.

Two sort paths:

- **`--low-mem-ssd`** (recommended for ≥ ~1 M reads): two-pass over a flat 16 B-per-read index (`lowMemNode`), with parallel chromosome-level processing in pass 2.
- **Classic** (default): loads input into a `seqread[]` array, sorts an index array with index-array `std::sort` (default) or bucket sort (opt-in via `--bucket-cutoff`).

### Input dispatch

All input paths land at the parser with `useMmap = true`:

- **File**: native `mmap` (PROT_READ | PROT_WRITE so the parser can NUL-terminate lines in place).
- **Stdin**: `fread`-slurped into a malloc'd buffer up front.
- **Gzip** (`.gz`): popened through `gzip -dc`, then slurped.

The parser then dispatches to either:

- `parseMmapDispatch` (classic path): newline-aligned chunks parsed concurrently into per-chunk per-chromosome partial linked lists, merged in a final serial step. `-t 1` short-circuits to `parseMmapSerial` (a thin wrapper around `parseLines<true>`) to skip the pre-count pass.
- `lowMemSortMmap` (low-mem path): similar parallel chunk parsing, but builds the 16 B-per-read `lowMemNode` table instead of `seqread[]`.

### Low-memory SSD path (`--low-mem-ssd`)

1. **Pass 1**: walk input, build `lowMemNode[]` (one 16 B entry per line: `size_t off; uint32_t next; int beg`) with per-chromosome linked lists. Parallel over chunks.
2. **Pass 2**: for each chromosome (in alphabetical order), walk its linked list, sort by `beg` (or `(beg, end)` for `--sort b`, or 5'-pos for `--sort 5`) with `std::sort`, then emit lines via `fwrite_unlocked` directly from mmap pointers. Pass 2 parallelises across chromosomes; output is serialised through a `mutex + condition_variable + nextChromToPrint` barrier so chromosomes appear in alphabetical order.

Per-chromosome output is a pre-sized `ChromBuf` (custom append-to-`char*` wrapper). v3.0.2 replaced the prior `open_memstream` to avoid its doubling-realloc churn at large chromosome sizes.

### Classic sort path

1. Walk per-chromosome linked lists, stamp `chrIdx` into the `seqread` union slot (overwriting `next`, no longer needed), collect indices into a flat `order[]` array. Track `chromStart[]` boundaries during the walk for the per-chrom radix path.
2. **Sort**:
   - `--sort s` and `--sort 5`: global LSD radix sort over packed 64-bit keys (high 32 = chrIdx, low 32 = position) when `n ≥ RADIX_SORT_THRESHOLD` (serial) or `RADIX_SORT_THRESHOLD_PAR` (parallel). Below threshold: `std::sort` with templated `ReadCmp<SortMode>` comparator (parallel exec policy at `-t > 1`).
   - `--sort b`: per-chromosome radix sort (`sortIndicesPerChromB`). Each chromosome's slice of `order[]` is sorted independently with a 64-bit `(beg, end)` key (the global `(chrIdx, beg, end)` would be 96 bits, too wide for `radixSort64`). Chromosomes processed in parallel via `std::for_each(par)` at `-t > 1`; serial radix per chrom.
3. **Emit** lines from `order[]`: full mmap line via `fputs_unlocked` (RAL or mmap path), or reconstruct `chr\tbeg\tend` + tail via `writeBedLine`.

`--bucket-cutoff` swaps step 2 for bucket sort. Body in `processChromBucketTpl<bool FullLine>` (`always_inline`, templated on full-line vs reconstruct-from-tail). Per chromosome, scatter reads into `chromTable[chosenPosition]`, then scan positions in order. Single-thread reuses one slab; multi-thread allocates per-thread per-chrom slabs and flushes through the producer-consumer barrier.

### Parallel parsing (mmap path, classic)

`parseMmapDispatch` splits the body into newline-aligned chunks (`parseChunkMmap` per chunk). Each chunk produces a `ChunkResult` with a per-chromosome partial linked list using **global** read indices, so no rebase pass is needed at merge time. The merge step concatenates per-chunk lists into the global `chrInfo`. `--collapse` allocates one `Arena` per chunk for weight strings; arenas are owned by `main` and freed at process exit.

### Constants

Top of `src/pioSortBed.cpp`:

- `kWeightBufSize` = 256 — stack buffer for BED score / RAL weight field. Over-long values rejected with a clear error.
- `kDefaultBucketChromBudget` = 4 GB — implicit per-chromosome `chromTable` rejection threshold. Override with `--max-mem`.
- `defaultBucketCutoff` = -1 — bucket sort disabled by default. Override with `--bucket-cutoff`.
- `RADIX_SORT_THRESHOLD` = 100 k, `RADIX_SORT_THRESHOLD_PAR` = 3 M — radix vs `std::sort` thresholds (in `sortIndices<>`).

### Limits (all dynamic / runtime / error-on-overflow)

The 3.0.x cleanup pass removed every previously hard-coded compile-time limit. Current state:

- **Line length**: no fixed limit. Stdin/gzip uses `getline` before slurp, post-slurp the parser uses `memchr`.
- **Chromosome name length**: no fixed limit. Stored as pointer+length into the line buffer.
- **Chromosome coordinate** (`int beg, int end`): signed 32-bit ⇒ ≤ 2.15 Gbp per single coordinate. Lift via `int64_t` widening.
- **Read count**: `uint32_t` indices in `seqread::next` / `lowMemNode::next` ⇒ ≤ 4.29 B reads on every sort path. Lift via `uint64_t` widening.
- **Per-chromosome bucket-sort RAM**: 4 GB default; `--max-mem N[GMK]` overrides.
- **BED score field length** (col 5): 255 bytes; over-long values error with a useful message.

### Repository layout

- `src/pioSortBed.cpp`, `src/CLI11.hpp` — source
- `test/test.sh` — 29-check correctness suite
- `benchmark/benchmark.sh` — cross-tool timing benchmark on synthetic BED6 fixtures (11 tools: `pio` at 1/4/8 threads, `pio-lm` at 1/4/8 threads, `sort` at 1/4/8 threads, `bedtools`, `bedops`)
- `benchmark/benchmark_na12878.sh` — real-data NA12878 BAM-derived benchmark
- `benchmark/bench_max_mem.sh` — `--max-mem` budget sweep on `pio-lm -t 4 @ 200 M`
- `benchmark/plot_readme.gp` — gnuplot driver for the four README plots (time/memory × log/linear)
- `benchmark/plot_maxmem.gp` — dual-y plot for the `--max-mem` sweep
- `benchmark/maxmem_sweep.csv`, `benchmark/maxmem_sweep.png` — sweep data + rendered plot
- `benchmark/history/` — per-version snapshots of `benchmark_readme.csv` (auto-archived by `benchmark.sh`); plots show the latest only, history kept for manual investigation
- `CHANGELOG.md` — version-keyed changelog (Keep-a-Changelog format)

## Memory profile

See the README's "Peak Memory (RSS)" section for measurements at 14 sizes (10 k – 200 M) across all configurations. Quick reference at 50 M reads (v3.0.8 numbers; the classic `-t N` rows were collected when bucket sort was the default at ≥ 50 M — from v3.0.10 onwards bucket sort is opt-in only):

- `-t 1` (classic, `std::sort`): ~4.1 GB peak — mmap input + `seqread[]`.
- `-t 8` (classic + `--bucket-cutoff 0`): ~12.6 GB peak — per-thread `chromTable` slabs sized to each chromosome's max coordinate.
- `--low-mem-ssd -t 1`: ~2.9 GB peak — lowest-RAM pioSortBed mode.
- `--low-mem-ssd -t 4`: ~3.9 GB peak (uncapped) — memory / speed sweet spot.
- `--low-mem-ssd -t 8`: ~4.8 GB peak (uncapped, default in `benchmark.sh` from v3.0.12). **Recommended fast path.**

`--max-mem` is a safety cap, not a memory optimiser. The mmap input + `lowMemNode` table dominate peak RSS; the flag only affects per-chrom `chromTable` / output buffer concurrent allocation. See README's `--max-mem` sweep section for the full curve.
