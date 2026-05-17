# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```bash
make            # produces ./pioSortBed (links libtbb.so dynamically)
make test       # runs the test suite under test/test.sh
                # (25 in the default build, 28 with WITH_HTSLIB=1,
                #  32 with WITH_HTSLIB=1 WITH_LOCISS=1)
make install    # installs to $PREFIX/bin (default /usr/local)
```

**Dependencies:** GCC ≥ 10 (C++20), oneTBB (`libtbb-dev` on Debian/Ubuntu).
CLI11 is bundled in `src/CLI11.hpp`.

The default build links libstdc++ and libgcc statically (`-static-libstdc++ -static-libgcc`) but leaves libtbb dynamic — Debian/Ubuntu only ships `libtbb.so`, not `libtbb.a`. For the GitHub release artefact (a fully self-contained binary) build a static oneTBB locally and use `make release-binary TBB_LIB=/path/to/libtbb.a`. The Makefile comments include the oneTBB cmake recipe.

Manual compilation:
```bash
g++ src/pioSortBed.cpp -Isrc -o pioSortBed -O3 -std=c++20 \
    -static-libstdc++ -static-libgcc -ltbb -llz4 -lzstd \
    -DVERSION_STRING=\"3.8.0\"
```

**Optional build flags** (all opt-in; the default `make` has zero
htslib/Arrow dependency):

```bash
make WITH_HTSLIB=1                                    # --bgzip, --tabix (uses system libhts)
make WITH_HTSLIB=1 HTSLIB=/path/to/htslib             # --bgzip, --tabix (statically linked)
make WITH_BAM=1 HTSLIB=/path/to/htslib                # BAM input (implies WITH_HTSLIB)
make WITH_RANS=1 HTSLIB=/path/to/htslib               # rans0/rans1 codecs for --external-merge (implies WITH_HTSLIB)
make WITH_LOCISS=1                                    # --lociss-output (Parquet)
make WITH_HTSLIB=1 WITH_LOCISS=1                      # both (independent)
```

- `WITH_HTSLIB=1` adds `--bgzip` and `--tabix` for integrated
  BGZF-compressed + indexed output. A drainer thread reads bytes off a
  stdout-redirected pipe and feeds them to `bgzf_write`; after BGZF
  close, `tbx_index_build3` builds the `.tbi` in place.
- `WITH_BAM=1` adds an in-RAM coordinate sort path that fires on
  `*.bam` input, reusing `radixSort64` and writing BAM via
  `sam_write1` with `hts_set_threads()` driving htslib's BGZF
  worker pool. Reads everything into memory (no spill); peak RSS
  is ~3× the BAM's decompressed size. Only the default coordinate
  sort is supported on BAM — `--collapse`, `--low-mem-ssd`,
  `--natural-sort`, and `--sort=b|5` are rejected.
- `WITH_LOCISS=1` adds `--lociss-output FILE` which writes a sorted
  Apache Parquet file per FORMAT_SPEC v2 (Chr, Start, End, typed
  per-flavor columns for BED4/5/6/12, MaxEndSoFar last; catch-all
  `Tail` string column for non-standard widths). Requires
  `libarrow-dev` + `libparquet-dev` (Debian/Ubuntu).

`HTSLIB` (when set) must point at a built htslib tree containing
`libhts.{a,so}` and the `htslib/` header directory.

## Usage

```bash
./pioSortBed file.bed > sorted.bed
./pioSortBed file.bed.gz > sorted.bed     # gzip transparently decompressed (slurped to memory)
./pioSortBed - < file.bed > sorted.bed    # stdin (also slurped)
./pioSortBed -o sorted.bed file.bed       # explicit output file
```

Options:
- `--sort s` (by start, default), `--sort b` (start + end), `--sort 5` (5'-end, strand-aware)
- `-n` / `--natural-sort` — `chr2 < chr10` instead of lexicographic
- `--collapse` — sum weights of regions sharing (chr, start)
- `--low-mem-ssd` — two-pass mode, **recommended fast path** for any input ≥ ~1 M reads
- `--external-merge` — streaming external sort for inputs > RAM (file input only). Pass 1 fills runs up to `--max-mem`, sorts each with `radixSort64`, writes binary compressed temp files; pass 2 k-way merges. Peak RSS bounded by `--max-mem`.
- `--multi-pass` — K-pass scan, **zero disk writes**. Pass 1 builds a `(chr, 1 MiB-quantum)` histogram + bin-packs into K groups ≤ `--max-mem`; passes 2..K+1 re-stream the input, filtering by group, sort + emit. Use for medium inputs (~1-5× RAM) where SSD-wear matters and you can afford `(K+1)×` reads. (file input only)
- `--merge-codec=raw|lz4|zstd|rans0|rans1` — temp-file codec for `--external-merge` (default zstd; rans0/rans1 require `WITH_RANS=1`)
- `--tmpdir=PATH` — temp-file location for `--external-merge` (default: `$TMPDIR` or `/tmp`)
- `--bgzip` (WITH_HTSLIB build only) — write BGZF-compressed BED text instead of plain text; requires `-o FILE`; mutually exclusive with `--lociss-output`
- `--tabix` (WITH_HTSLIB build only) — after `--bgzip`, build a `.tbi` index in place; implies `--bgzip`
- `--lociss-output FILE` (WITH_LOCISS build only) — write a sorted Apache Parquet file per FORMAT_SPEC v2 instead of BED text. Pairs with `--collapse` on the classic sort path (writes a five-column `{Chr, Start, End, Score double, MaxEndSoFar}` schema per FORMAT_SPEC §10); still rejects `--sort=b|5`
- `-t N` / `--threads N` — worker pool size (`0` = all cores, default; `1` = serial)
- `--max-mem=N[GMK]` — concurrent per-chromosome emit-buffer cap on the `--low-mem-ssd` parallel pass 2, OR per-run memory budget on `--external-merge` (default 1 GiB on the latter)
- `-o FILE` / `--output FILE` — write to file (default: stdout)
- `-v` / `--verbose` — opt in to timing / chromosome-length info on stderr (default: silent)

BED header lines (`#`, `track `, `browser `) pass through unchanged.

## Architecture

Single source file (`src/pioSortBed.cpp`, ~5600 lines). Section banners divide it into:
**CORE TYPES & ALLOCATOR** / **PARSERS** / **LOW-MEMORY SSD PATH** / **CLASSIC SORT PATH** / **Parallel mmap parsing** / **BAM PATH** / **EXTERNAL MERGE SORT PATH** / **MULTI-PASS NO-WRITES PATH** / **LOCISSD OUTPUT** / **CLI / MAIN**. Grep `^// =\+$` to jump between sections.

Four sort paths:

- **`--low-mem-ssd`** (recommended for inputs that fit in RAM, ≥ ~1 M reads): two-pass over a flat 16 B-per-read index (`lowMemNode`), with parallel chromosome-level processing in pass 2.
- **Classic** (default): loads input into a `seqread[]` array, sorts an index array with `radixSort64` over packed `(chrIdx, pos)` 64-bit keys, falls back to `std::sort` below a small-input threshold.
- **`--multi-pass`** (zero writes, ~1-5× RAM regime): pass 1 streams via `getline`, builds a `(chr, beg-quantum)` histogram, bin-packs into K consecutive groups ≤ `--max-mem`. Passes 2..K+1 re-stream filtering by group, sort + emit. Total reads = (K+1)×file, total writes = 0. SSD-wear-friendly fallback when you don't want temp files.
- **`--external-merge`** (for inputs > RAM): pass 1 streams via `getline`, fills a per-run buffer up to `--max-mem`, sorts with `radixSort64`, writes a binary compressed temp file (codec via `--merge-codec`, default zstd). Pass 2 is a k-way min-heap merge over all run files. Peak RSS bounded by `--max-mem` regardless of input size. zstd is fastest at scale even vs uncompressed temps because the I/O saved exceeds the compression CPU cost; rans0 (htscodecs, `WITH_RANS=1` only) is even faster on this workload thanks to SIMD-vectorised rANS.

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
3. **Emit** lines from `order[]`: full mmap line via `fputs_unlocked` (mmap path), or reconstruct `chr\tbeg\tend` + tail via `writeBedLine` (stdin/gzip path, where only the tail is in the arena).

### Parallel parsing (mmap path, classic)

`parseMmapDispatch` splits the body into newline-aligned chunks (`parseChunkMmap` per chunk). Each chunk produces a `ChunkResult` with a per-chromosome partial linked list using **global** read indices, so no rebase pass is needed at merge time. The merge step concatenates per-chunk lists into the global `chrInfo`. `--collapse` allocates one `Arena` per chunk for weight strings; arenas are owned by `main` and freed at process exit.

### Constants

Top of `src/pioSortBed.cpp`:

- `kWeightBufSize` = 256 — stack buffer for the BED weight field on the `--collapse` path. Over-long values rejected with a clear error.
- `RADIX_SORT_THRESHOLD` = 100 k, `RADIX_SORT_THRESHOLD_PAR` = 3 M — radix vs `std::sort` thresholds (in `sortIndices<>`).

### Limits (all dynamic / runtime / error-on-overflow)

The 3.0.x cleanup pass removed every previously hard-coded compile-time limit. Current state:

- **Line length**: no fixed limit. Stdin/gzip uses `getline` before slurp, post-slurp the parser uses `memchr`.
- **Chromosome name length**: no fixed limit. Stored as pointer+length into the line buffer.
- **Chromosome coordinate** (`int beg, int end`): signed 32-bit ⇒ ≤ 2.15 Gbp per single coordinate. Lift via `int64_t` widening.
- **Read count**: `uint32_t` indices in `seqread::next` / `lowMemNode::next` ⇒ ≤ 4.29 B reads on every sort path. Lift via `uint64_t` widening.
- **BED score field length** (col 5): 255 bytes; over-long values error with a useful message.

### Repository layout

- `src/pioSortBed.cpp`, `src/CLI11.hpp` — source
- `test/test.sh` — correctness suite (25 / 28 / 32 across build configurations)
- `.github/workflows/ci.yml` — 3-config matrix CI on every push and PR
- `benchmark/benchmark.sh` — cross-tool timing benchmark on synthetic BED6 fixtures (11 tools: `pio` at 1/4/8 threads, `pio-lm` at 1/4/8 threads, `sort` at 1/4/8 threads, `bedtools`, `bedops`)
- `benchmark/benchmark_na12878.sh` — real-data NA12878 BAM-derived benchmark
- `benchmark/bench_max_mem.sh` — `--max-mem` budget sweep on `pio-lm -t 4 @ 200 M`
- `benchmark/bench_query.sh` + `bench_query_drivers.py` — region-query latency benchmark (4 readers × 3 sizes × 2 cache states)
- `benchmark/bench_pipeline.sh` — end-to-end pipeline wall time + bytes written
- `benchmark/bench_loci.sh` — Loci (Python; polars + pandas) sort scaling sweep across the same sizes as `benchmark.sh`; needs `LOCI_PY`/`LOCI_SORT_SCRIPT` (defaults match the bench desktop layout)
- `benchmark/plot_readme.gp` — gnuplot driver for the four README plots (time/memory × log/linear)
- `benchmark/plot_maxmem.gp` — dual-y plot for the `--max-mem` sweep
- `benchmark/maxmem_sweep.csv`, `benchmark/maxmem_sweep.png` — sweep data + rendered plot
- `benchmark/history/` — per-version snapshots of `benchmark_readme.csv` (auto-archived by `benchmark.sh`); plots show the latest only, history kept for manual investigation
- `CHANGELOG.md` — version-keyed changelog (Keep-a-Changelog format)
- `paper/preprint.md`, `paper/references.bib`, `paper/figures/` — preprint source

## Memory profile

See the README's "Peak Memory (RSS)" section for measurements at 14 sizes (10 k – 200 M) across all configurations. Quick reference at 50 M reads:

- `-t 1` (classic, `radixSort64`): ~4.1 GB peak — mmap input + `seqread[]`.
- `--low-mem-ssd -t 1`: ~2.9 GB peak — lowest-RAM pioSortBed mode.
- `--low-mem-ssd -t 4`: ~3.9 GB peak (uncapped) — memory / speed sweet spot.
- `--low-mem-ssd -t 8`: ~4.8 GB peak (uncapped, default in `benchmark.sh`). **Recommended fast path.**

`--max-mem` is a safety cap on the `--low-mem-ssd` parallel pass-2 emit-buffer budget, not a memory optimiser. The mmap input + `lowMemNode` table dominate peak RSS regardless. See README's `--max-mem` sweep section for the full curve.
