# Changelog

All notable changes to pioSortBed are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project uses [semantic versioning](https://semver.org).

## [3.0.2] — 2026-04-30

### Changed
- **`--low-mem-ssd` pass 2 now writes through a pre-sized `ChromBuf`** instead
  of a glibc `open_memstream`. The memstream's exponential-doubling realloc
  was holding both old and new buffers transiently at every growth step,
  so peak RSS was inflating ~2× the chromosome's final output size during
  the doubling. Pre-sizing to `chromCount × avgLineBytes × 1.05` eliminates
  the doublings on the parallel path.
  - Single-thread (`-t 1`): keeps a 64 KB scratch reused across chromosomes,
    `fwrite_unlocked`'d to stdout when full. Peak RAM unchanged; ~3% faster.
  - Parallel (`-t > 1`): pre-sized per-chromosome buffer, queued through the
    existing producer-consumer barrier. **20M reads measured: 2.32 s / 2.31 GB
    → 1.84 s / 1.85 GB — 21% faster, 20% less RAM.**
- `processChrom` now takes a `ChromBuf&` instead of a `FILE*`; collapse mode
  uses `vsnprintf` into a 512-B stack buffer + `memcpy` instead of `fprintf`.

## [3.0.1] — 2026-04-30

### Changed
- **`--low-mem-ssd` node table shrunk from 24 B to 16 B per read** by dropping
  the cached `end` field. The default `--sort s` and `--collapse` modes never
  needed it; `--sort 5` already re-parses for strand and pulls `end` out of
  that same call for free; `--sort b` now re-parses `end` from the line
  (one tab-scan + atoi per read in pass 2 — cheap).
  - 50M reads: ~0.4 GB saved on the node table (1.2 GB → 0.8 GB)
  - 200M reads: ~1.6 GB saved (4.8 GB → 3.2 GB)
  - As a bonus, the smaller cache footprint speeds pass 2 up: 20M-read
    `--low-mem-ssd -t 1` measured at 5.46 s / 1.18 GB (was 6.39 s / 1.30 GB
    in v3.0.0) — ~15% faster, ~9% lower peak RSS at the same input.

---

## [2.2.0] — 2026-04-29

### Added
- **`--max-mem=N[GMK]`** — memory budget for the parallel bucket-sort path.
  Caps concurrent per-chromosome `chromTable` allocations via a
  `mutex+condition_variable` gate; tasks acquire their `chromTable`'s
  byte cost from a shared budget before starting and release after the
  chromosome is fully processed. A chromosome whose own cost exceeds the
  budget runs alone (clamping its draw to the full budget). Default:
  no cap (preserves the v2.1.x behavior).

### Performance
- v2.1.1: LSD radix sort for `--sort s/5` at -t 1 (~2× at 10M).
- v2.1.2: parallel LSD radix sort for `--sort s/5` at -t > 1 (~8% at 5M+).
- v2.1.3 / 2.2.0: optional `--max-mem` budget cap (see Added above).

### Other improvements since v2.1.0
- 752ae9c: -t 1 `fullLine` hoist (recovers a small post-integration regression).
- 29d1131: `MAP_POPULATE` + `MADV_HUGEPAGE` on the mmap path (cold-cache win).
- daba255: `benchmark/history/` directory for per-version CSV snapshots; benchmark
  SIZES extended to include 200k, 2M, 20M.
- bad5b2c: `unordered_map<string, V>` → flat-vector `ChrNameMap<V>` template
  (-7% instructions, -8.6% branches, wall time within noise).

### Memory cap measurements (50M reads, 10-chrom benchmark fixture, -t 8)
| `--max-mem` | wall time | peak RSS | vs no cap |
|------------:|----------:|---------:|----------:|
| (no cap)    | 6.37 s    | 12.6 GB  | baseline  |
| `8G`        | 5.88 s    | 12.6 GB  | tied      |
| `6G`        | 8.16 s    | 10.2 GB  | -19% RAM, +28% time |
| **`4G`**    | **8.17 s** | **8.6 GB** | **-32% RAM, +28% time** |
| `2G`        | 12.06 s   | 6.5 GB   | -48% RAM, +89% time |

Sweet spot is roughly `--max-mem=4G` on this hardware: a 32% peak-RAM cut
for a 28% wall-time hit. Pick a budget that fits your RAM minus
~3 GB for the input mmap, reads array, and runtime overhead.

---

## [2.2.1] — 2026-04-29

### Performance
- **`--collapse` no longer disables parallel mmap parsing.** Previously, when
  the input was an mmap'd file with `--collapse` and `--threads > 1`,
  `parseMmapDispatch` fell back to the legacy single-threaded
  `parseMmapSerial` because the per-read weight strings need to be copied
  out of the mmap (the line buffer gets `\n` overwritten with `\0` during
  parse). The parallel parser now allocates one `Arena` per chunk so each
  thread copies its weight strings into its own buffer with no
  cross-thread contention, and the collapse path runs at full parallelism.

  Measured on 10M random BEDWEIGHT rows:
  - `-t 1`: 5.72 s (`Reading has taken 1 seconds`)
  - `-t 4`: 5.15 s (`Reading has taken 0 seconds` — parsing accelerated)
  - `-t 8`: 5.14 s (same)

  ~10% wall-time win at 10M; bigger win expected at 50M+ where parsing
  dominates more. Output phase is unchanged (still single-threaded
  collapsed-output writer).

### Notes
- Per-chunk arenas are owned by `main` (a `std::vector<Arena*>` passed
  to `parseMmapDispatch` by reference) so they outlive the parser and
  remain valid through the bucket-sort/output phase, where
  `reads[i].line` points into them.

---

## [2.2.2] — 2026-04-29

### Performance
- **`--low-mem-ssd`: ~26% faster, ~3.5% lower peak RSS at 50M reads.**
  Two changes, made together so the wins stack:
  1. **`uint32_t` line offsets in `lowMemNode`.** Previously `size_t off`
     used 8 bytes; the `next` field's 4-byte alignment forced 4 bytes of
     padding for a 16-byte node. Switching to `uint32_t off` opens the
     freed 8 bytes for inline `beg`+`end` storage at no memory cost (still
     16 B/node). Adds a sanity check that bails out cleanly if the input
     mmap exceeds 4 GB (BED equivalent: well above 100M reads).
  2. **Pre-parsed `beg`+`end` in `lowMemNode`.** Pass 2 used to re-parse
     every line a second time to populate `lowMemRec`. With `beg`+`end`
     already in the node, `--sort s` and `--sort b` skip pass-2 parsing
     entirely and sort a contiguous 12-byte `{beg, end, off}` array (vs
     the prior 24-byte `lowMemRec` — also halves the per-chromosome sort
     buffer). `--sort 5` and `--collapse` still fetch one extra field
     (strand or weight respectively) but no longer waste a parse on
     `chr`/`beg`/`end`.

  Measured (50M reads, 10-chrom benchmark fixture, perf stat -r 5):

  |                | before | v2.2.2 |
  |----------------|-------:|-------:|
  | task-clock     | 18.8 s | 13.5 s |
  | cycles         | 77.6 B | 53.6 B |
  | instructions   | 71.8 B | 53.1 B |
  | wall time      | 10.2 s | **7.54 s** (-26%) |
  | peak RSS       | 3.15 GB | **3.04 GB** (-3.5%) |

### Notes
- Output is bit-identical to the previous `--low-mem-ssd` (verified by
  diffing both runs' stdout on a 50M-row fixture).

---

## [2.2.3] — 2026-04-29

### Performance
- **`--low-mem-ssd` pass 1 now parses in parallel.** Previously pass 1 was
  always single-threaded; pass 2's per-chrom sort already used parallel
  `std::sort`. Pass 1 dominates total time on cold inputs, so going parallel
  there gives the biggest remaining `--low-mem-ssd` lever.

  The parallel chunk parser mirrors the regular `parseMmapDispatch`: split
  the body of the mmap into newline-aligned chunks, count newlines per
  chunk to size `nodes[]` exactly, parse each chunk into its slice using
  global node indices (per-chunk ChrNameMap partials), then merge by
  prepending each chunk's per-chrom list to the global one. Per-chunk
  partials are O(K_chroms) total, ~80 KB at 22 chunks — the per-thread
  memory cost is genuinely negligible, unlike the regular path's per-chrom
  `chromTable` slabs.

  Header pre-pass extracted: leading `#`/`track`/`browser` lines are now
  emitted once before chunking, matching the regular path's behaviour.
  Mid-file headers are no longer passed through (they would have failed
  to parse as data anyway).

  Measured (50M reads, perf stat -r 5):

  |     | wall time | peak RSS |
  |----:|----------:|---------:|
  | -t 1 | 10.05 s   | 2.99 GB  |
  | -t 4 | 6.58 s    | 3.03 GB  |
  | -t 8 | **5.91 s** | **3.03 GB** |

  At `-t 8`, `--low-mem-ssd` is now essentially as fast as the regular
  parallel path (6.15 s at the same input) but at **~4× lower peak RAM**
  (3.0 GB vs 12.6 GB). The previous "low-mem mode is for the memory-tight
  scenario, accept the speed cost" trade-off is mostly gone.

### Notes
- Output line order between `-t 1` and `-t 8` may differ on tie-broken
  reads (sort is unstable across thread counts), but the multiset is
  identical and the sort key remains monotonic — same correctness
  semantics as the regular `-t > 1` paths.

---

## [2.2.4] — 2026-04-29

### Performance
- **`--low-mem-ssd` pass 2 now runs in parallel across chromosomes.** With
  parallel pass 1 already in v2.2.3, this completes the parallelisation:
  per-chromosome work runs through `std::for_each(std::execution::par, ...)`,
  each chromosome writes to an `open_memstream` buffer, and a producer-
  consumer barrier flushes buffers to stdout in alphabetical order.
  `--max-mem=N[GMK]` caps concurrent per-chromosome buffers (rough cost
  estimate: ~62 B per read = 12 B sort-rec + ~50 B output buffer).

### `--low-mem-ssd` is now the recommended path

50M reads, -t 8, 10-chrom benchmark fixture (median of 3):

| Mode | Wall time | Peak RSS |
|------|----------:|---------:|
| regular `-t 8` | 6.47 s | 12.6 GB |
| `--low-mem-ssd -t 8` | **3.63 s** | **5.4 GB** |
| `--low-mem-ssd -t 8 --max-mem=2G` | 4.35 s | 5.0 GB |

`--low-mem-ssd` is now **1.78× faster AND 2.3× lower peak RSS** than the
regular path. Strict win on both axes for mmap input. Per-thread state
is bounded by output-buffer cost (~250 MB max for chr1) instead of the
regular path's per-chromosome `chromTable` slabs (up to 1 GB each).

The regular path remains the fallback for stdin/gzip input, where mmap
isn't available. For everyday file-input usage you almost always want
`--low-mem-ssd` now.

### Notes
- Output line order between `-t 1` and `-t 8` may differ on tie-broken
  reads (parallel sort is unstable across thread counts), but the
  multiset and sort key are identical — same semantics as the regular
  `-t > 1` paths.
- `--max-mem=4G` is essentially no-op on this fixture because the
  natural peak is already 5.4 GB; pick a tighter budget to actually
  trade speed for memory.

---

## [3.0.0] — 2026-04-30

Major version bump signalling the cumulative shift since 2.0.x: `pioSortBed`
is now markedly faster than every other tool we tested across the full size
range, and `--low-mem-ssd` has gone from "use this when RAM is tight" to
"use this by default" — it's now strictly the fastest configuration for any
mmap-able file >50k reads.

### Bug fix
- **`--low-mem-ssd` no longer caps the input file at <4 GB.** v2.2.2's
  `uint32_t off` packing in `lowMemNode` (intended to be memory-neutral)
  silently bailed out at files ≥4 GB, which silently aborted the benchmark
  when it reached 100M reads. Reverted to `size_t off`; node size grows
  from 16 B to 24 B (the trailing `int` requires 4 B padding for `size_t`
  alignment), which is fine: the speed wins from pre-parsed `beg`+`end`
  (also v2.2.2) are unaffected, and at 200M reads the node table is
  4.8 GB — well within the 21 GB peak.

### New benchmark coverage (100M and 200M reads)
The benchmark now runs through 200M-row BED6 fixtures (~8.6 GB input).
At those sizes `pioSortBed -t 8` and `bedtools` are skipped because
their memory growth (per-thread `chromTable` slabs and linear-in-input
respectively) would exceed the 30 GB RAM available on the test box.

#### 200M reads — wall time
| Tool | Wall | Peak RSS |
|---|---:|---:|
| **`pioSortBed --low-mem-ssd`** (default `-t 0`, `--max-mem=4G`) | **22.72 s** | 21.1 GB |
| GNU sort `--parallel=8` | 1min 45.5 s | 15.4 GB |
| bedops `sort-bed` | 2min 34.6 s | 10.4 GB |
| `pioSortBed -t 1` | 1min 8.4 s | 13.6 GB |
| GNU sort `--parallel=1` | 5min 28.5 s | 15.4 GB |

`pioSortBed --low-mem-ssd` is **4.6× faster than GNU sort 8t**, **6.8×
faster than bedops**, and **14.5× faster than GNU sort 1t** at 200M.

### Tooling
- `benchmark.sh` extended to 14 size points (10k → 200M).
- For huge sizes (`>= 100M`):
  - `pio8` and `bedtools` skipped (would OOM the 30 GB system).
  - `pio-lm` runs with `--max-mem=4G` to cap concurrent per-chromosome buffers.
  - GNU sort uses `--buffer-size=50%` (down from 80%) so its in-process
    buffer doesn't race against `pio-lm`'s working set.
  - Tool outputs go to `/dev/null` instead of the tmpfs `$TMPDIR/output.tmp`
    so multi-GB outputs don't fill the workdir.
- Recommend running with `TMPDIR=/var/tmp bash benchmark/benchmark.sh` for
  huge fixtures: GNU sort's external-sort spillover and the 8 GB BED file
  together exceed the typical 16 GB tmpfs ceiling.

### Notes
- `pio-lm` at 100M+ exclusively uses `--max-mem=4G` in the benchmark; without
  that cap the default-thread-count parallelism would queue ~20 GB of
  per-chromosome output buffers.
- The 200M `pio-lm` measurement is technically run with default-but-capped
  parallelism. Without the cap, peak RSS at -t 0 (22 cores) on this fixture
  would push past 30 GB and risk OOM.

> **Post-3.0.0 doc update (2026-04-30):** the `pio-lm` row above was
> measured with the default thread count (= all cores, 22 on the bench box),
> while the rest of the table compared `-t 1` and `-t 8` configs. The
> README and `benchmark/benchmark_readme.csv` were since re-measured with
> explicit `-t 1` and `-t 8` low-mem variants for an apples-to-apples
> comparison. The 22-thread `pio-lm` numbers above are kept here as
> historical record; the README values supersede them. Headlines after
> the re-measurement: `pio-lm -t 8` at 200M is **35.16 s** (3.0× GNU sort
> 8t, 4.4× bedops, 9.3× GNU sort 1t).

---

## [Unreleased]

### Performance
- `processChromBucket` re-templated on `bool FullLine` so the per-read
  inner loops in the bucket-sort path are fully specialized at compile
  time (gcc couldn't constant-fold the runtime parameter through the
  function call boundary). Hot-cache 50M / `-t 1` median wall time on a
  6-chromosome fixture: 14.18 s baseline → 13.36 s (~6% speedup, recovers
  the small regression vs origin's pre-integration baseline).
- mmap path now passes `MAP_POPULATE` (pre-fault all pages at mmap time)
  and `madvise(MADV_HUGEPAGE)` (hint transparent huge pages) — both Linux
  extensions, gated on `#ifdef`. Mostly a cold-cache improvement (large
  files no longer fault-in millions of pages incrementally during the
  parse pass).

### Tooling
- `benchmark/benchmark.sh` extended with three more data points (200k, 2M,
  20M) for finer resolution. After each run it archives a versioned
  snapshot of `benchmark_readme.csv` to `benchmark/history/` so the per-
  release numbers are preserved (only the latest is plotted).

### Benchmark refresh (numbers in README updated)
50M / `-t 8` is now **9.27 s** vs `25.94 s` for GNU sort 8t (3.1×),
`38.50 s` for bedops (4.1×), and `1min06.6s` for bedtools (7.2×). New
mid-range data points reveal `--low-mem-ssd` is actually the fastest
pioSortBed mode at 10M–20M (the parallel-sort scaffolding overhead
doesn't fully pay off in that range).

---

## [2.1.0] — 2026-04-28

### Added — features
- **BED header passthrough.** Lines starting with `#`, `track `, or `browser `
  are emitted to stdout unchanged (BED convention places these at the top of
  the file, before any sortable rows).
- **Gzip-compressed input.** Files with `.gz` extension are transparently
  decompressed via `gzip -dc | popen()`. The stdin/gzip parser path uses the
  arena allocator; mmap path is unchanged.
- **`--natural-sort` / `-n`** — natural chromosome ordering (`chr2 < chr10`
  instead of lexicographic `chr10 < chr2`). Embedded digit runs are compared
  numerically. Default behaviour (lexicographic) is unchanged.

### Added — performance
- **Parallel mmap parsing.** Input is split into newline-aligned chunks
  parsed concurrently into per-chunk per-chromosome partial linked lists,
  merged serially at the end. `-t 1` short-circuits to the legacy single-pass
  parser to avoid the pre-count overhead. Header lines are emitted to stdout
  in a leading pre-pass before chunking.
- **Parallel per-chromosome bucket sort.** Chromosomes processed
  concurrently via `std::for_each(std::execution::par, ...)`. Each thread
  allocates its own `chromTable[thisChrLen+1]` sized exactly to that
  chromosome. Per-chromosome output is buffered through `open_memstream` and
  flushed under a producer-consumer barrier so stdout sees chromosomes in
  alphabetical order while live buffer count stays bounded by the worker pool.
- **Test suite** in `test/test.sh` — 29 checks covering BED3/BED6/RAL,
  classic & forced-bucket sort, `--sort b/5`, `--collapse`, `--low-mem-ssd`,
  file/stdin/gzip paths, header passthrough, natural-sort, parallel paths
  (`-t 8`), and multi-chromosome stress. `make test` runs the suite.

### Changed
- **Parallelism backend: OpenMP / `__gnu_parallel::sort` → C++17
  `std::execution::par` over oneTBB.** `tbb::global_control` caps the worker
  pool to `--threads N`.
- **Build flags**: `-std=c++11 -fopenmp -lgomp` →
  `-std=c++17 -static-libstdc++ -static-libgcc -ltbb`. C++/GCC runtimes are
  statically linked; libtbb stays dynamic (Debian ships no `libtbb.a`).
- **`--threads`** option now also governs parallel parsing and parallel
  bucket sort (previously only the classic-sort path).
- **Output buffer**: 1 MB → 8 MB. Amortises write syscalls on 100M+ read
  files.
- **Refactors**:
  - `ReadCmp<SortMode>` templated comparator collapses six duplicated
    `std::sort` callsites in the classic-sort path.
  - `printRead` helper eliminates duplicated stdout-vs-rebuild branching.
  - `sumList` / `sumList2` merged into `parseWeight` + `sumWeightsBuf` +
    `sumWeightsList`; `strtof` replaces `sscanf` for weight parsing
    (faster + locale-independent).
  - `processChromBucket` extracted as a free function used by both serial
    and parallel bucket-sort paths.
  - `madvise(MADV_RANDOM)` before low-mem-ssd pass 2 (random access pattern).
- **Repository layout**: source moved to `src/`, benchmark to `benchmark/`,
  tests to `test/`. `Makefile` adds `make install` (`PREFIX=/usr/local`),
  `make test`, and `VERSION` macro injection via `-DVERSION_STRING`.

### Fixed
- **`--sort 5` was sorting by `beg` instead of 5'-end.** `parseBedLineFull`
  was double-skipping the delimiter after the weight field, so strand was
  parsed as NUL on every line and the strand-aware 5'-end calculation
  collapsed to plain `beg`.
- `realloc` null check in `parseLines` (OOM now prints an error instead of
  segfaulting on a NULL pointer).
- Comment `// 500Mb` updated to reflect the actual 1 Gbp limit.

### Performance — 50M-row BED6 fixture, warm cache, median of 3 runs
| Threads | v2.0.x | v2.1.0 | Speedup |
|--------:|-------:|-------:|--------:|
| `-t 1`  | 14.0 s | 14.7 s | 0.95×   |
| `-t 4`  | 13.9 s | 7.4 s  | **1.88×** |
| `-t 8`  | 14.0 s | 4.2 s  | **3.33×** |

Cross-tool at 50M / `-t 8`: pioSortBed is **3.2×** faster than GNU sort
`--parallel=8`, **5.5×** faster than `bedops sort-bed`, and **8.6×** faster
than `bedtools sort` on the same machine. See README for the full grid
(100k → 50M).

### Memory trade-off
Peak RSS at 50M / `-t 8` rises ~4 GB → ~13 GB because each thread allocates
its own per-chromosome `chromTable` and queued per-chromosome output buffers
accumulate. Use `-t 1` if RAM is tight; the single-threaded path keeps a flat
~4 GB ceiling.

### Notes
- `--collapse` on mmap input falls back to the serial parser at `-t > 1`
  (parallel parsing would need per-thread arenas for the weight strings).
- Cross-tool benchmark tables in the README and `benchmark_readme.csv` are
  refreshed for v2.1.0; older v2.0.0 numbers are preserved in git history.

---

## [2.0.1] — 2026-03-25
- Updated benchmark numbers and plots to current binary.

## [2.0.0]
- Added low-memory SSD mode (`--low-mem-ssd`): two-pass algorithm that keeps
  only line offsets in RAM and processes one chromosome at a time.
- Integrated benchmarking with cross-tool comparison.

## [1.2.1]
- Restored fast `fputs` output on the mmap path; added `--version` flag.

## [1.2.0]
- Parsing optimizations; tail-only line storage for stdin BED to reduce
  arena memory.

## [1.1.0]
- Benchmark results and plot generation with real data.

## [1.0.0]
- Major performance optimizations and modernization (CLI11, mmap input,
  hybrid classic/bucket sort strategy with `--bucket-cutoff`).
