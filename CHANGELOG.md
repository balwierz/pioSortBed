# Changelog

All notable changes to pioSortBed are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
the project uses [semantic versioning](https://semver.org).

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
