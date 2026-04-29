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
    -static-libstdc++ -static-libgcc -ltbb -DVERSION_STRING=\"2.2.1\"
```

## Usage

```bash
./pioSortBed file.bed > sorted.bed
./pioSortBed file.bed.gz > sorted.bed     # gzip input transparently decompressed
./pioSortBed - < file.bed > sorted.bed    # stdin
```

Options:
- `--sort s` (by start, default), `--sort b` (start + end), `--sort 5` (5'-end, strand-aware)
- `-r` / `--ral` — RAL input format instead of BED
- `-n` / `--natural-sort` — `chr2 < chr10` instead of lexicographic
- `--collapse` — sum weights of regions sharing (chr, start)
- `--low-mem-ssd` — two-pass file-input mode, peak RAM ≈ largest chromosome
- `--bucket-cutoff N` — switch to bucket sort at N reads (default 50M; 0 = always bucket)
- `-t N` / `--threads N` — worker pool size (0 = all cores; 1 = serial)
- `--max-mem=N[GMK]` — cap concurrent per-chromosome `chromTable` allocations in the parallel bucket-sort path

BED header lines (`#`, `track `, `browser `) pass through unchanged.

## Architecture

Single source file (`src/pioSortBed.cpp`, ~2000 lines). Hybrid sort: files with fewer than `--bucket-cutoff` reads (default 50M) take the classic O(n log n) path; larger files take the bucket/counting sort path with O(n + m) complexity (m = max chromosome length).

### Data flow
1. **Input dispatch.** mmap'd file → `parseMmapDispatch`. Stdin or `.gz` → arena-backed `parseLines<false>`.
2. **Header pre-pass** (mmap only). Leading `#`/`track `/`browser ` lines are emitted to stdout before chunking.
3. **Parsing** stores reads into a flat `seqread[]` array; per-chromosome linked lists thread through `seqread.next`. The `chrInfoT.lastRead` head is kept in a small flat-vector `ChrNameMap<chrInfoT>` (~70-line replacement for `unordered_map<string, V>`, faster on the miss path because it avoids `std::string` construction).
4. **Sort dispatch** chooses classic vs bucket sort based on `totalReads` vs `--bucket-cutoff`.
5. **Output.** Plain mode: `printRead` streams either the full mmap line (`fputs_unlocked`) or reconstructs from the stored "tail" (`writeBedLine`). Collapse mode: one `fprintf` per unique (chr, position) summing weights.

### Parallel parsing (mmap path)
`parseMmapDispatch` splits the body of the mmap into newline-aligned chunks and parses each in parallel via `std::for_each(std::execution::par, ...)`. Each chunk produces a `ChunkResult` with a per-chromosome partial linked list (using *global* read indices, so no rebase pass at merge time). The merge step concatenates per-chunk lists into the global `chrInfo`. `--collapse` allocates one `Arena` per chunk for weight strings; arenas are owned by `main` and freed at process exit.

`-t 1` short-circuits to `parseMmapSerial` (a thin wrapper around the legacy `parseLines<true>` realloc-grow parser) to skip the pre-count pass.

### Classic sort path
1. Walk per-chromosome linked lists, stamp `chrIdx` into the union slot (overwriting `next`, no longer needed), collect indices into a flat `order[]` array.
2. **Sort by `(chrIdx, position)`**:
   - `--sort s` and `--sort 5` with n ≥ 100k (serial) or n ≥ 3M (parallel): **LSD radix sort** over packed 64-bit keys (high 32 = chrIdx, low 32 = position). 8-bit digits, constant-byte detection skips passes where all keys agree on a byte. Parallel variant: thread-local histograms in parallel, serial 2D prefix sum over (T × 256), parallel scatter into disjoint per-(thread, bucket) ranges.
   - `--sort b` (any size) and small inputs: `std::sort` with templated `ReadCmp<SortMode>` comparator (`if constexpr` collapses each instantiation to a single int compare per arm at compile time). `-t > 1` uses `std::sort(std::execution::par, ...)`.
3. Linear scan to print, with `fullLine` hoisted out of the per-read loop (avoids a per-read branch).

### Bucket sort path
Per chromosome, scatter reads into `chromTable[chosenPosition]` (a sparse array indexed by position), then scan positions in order to emit reads. Body factored into `processChromBucketTpl<bool FullLine>` (templated on full-line vs reconstruct-from-tail; `__attribute__((always_inline))` so `gcc` constant-folds the FullLine branch into the inner loop).

- **Single-thread** (`-t 1`): one shared `chromTable[maxChrLen+1]` is reused across chromosomes; slots zero themselves during the consume pass.
- **Multi-thread** (`-t > 1`): chromosomes processed concurrently via `std::for_each(par, ...)`. Each thread allocates its own `chromTable[thisChrLen+1]` sized exactly to its chromosome. Per-chromosome output is written to an `open_memstream` buffer and flushed under a producer-consumer barrier (`mutex + condition_variable + nextChromToPrint`) so stdout sees alphabetical order while live buffer count stays bounded by the worker pool.
- `--max-mem=N[GMK]` adds a second `mutex+condvar` gate that caps concurrent `chromTable` allocations to the byte budget. Default: no cap (preserves previous behavior).

### Key compile-time constants (top of `src/pioSortBed.cpp`)
- `lineBufSize` = 1024 — max BED line length (stdin only; no limit for mmap)
- `chrNameBufSize` = 256 — max chromosome name length
- `chrLenLimit` = 1,000,000,000 — max coordinate value (1 Gbp); increase + recompile for non-standard genomes
- `defaultBucketCutoff` = 50,000,000 — hybrid sort threshold (overridden by `--bucket-cutoff`)
- `RADIX_SORT_THRESHOLD` = 100,000 — single-thread radix vs `std::sort`
- `RADIX_SORT_THRESHOLD_PAR` = 3,000,000 — multi-thread radix vs `std::sort(par)`

### Repository layout
- `src/pioSortBed.cpp`, `src/CLI11.hpp` — source
- `test/test.sh` — 29-check correctness suite
- `benchmark/benchmark.sh` — cross-tool timing benchmark on synthetic BED6 fixtures
- `benchmark/benchmark_na12878.sh` — real-data NA12878 BAM-derived benchmark
- `benchmark/plot_readme.gp` — gnuplot driver for the four README plots (time/memory × log/linear)
- `benchmark/history/` — per-version snapshots of `benchmark_readme.csv` (auto-archived by `benchmark.sh`); plots only show the latest, history is preserved for manual investigation
- `CHANGELOG.md` — version-keyed changelog (Keep-a-Changelog format)

## Memory profile

See the README's "Peak Memory (RSS)" section for measurements at 12 sizes (10k–50M) across all configurations. Quick reference at 50M reads:

- `-t 1` (any sort mode): ~4.1 GB peak (mostly mmap input + reads array; classic sort doesn't allocate per-chromosome scratch).
- `-t 8` (bucket sort): ~12.6 GB peak — per-thread `chromTable` slabs sized to each chromosome's max coord, plus queued per-chromosome output buffers. Use `--max-mem=4G` to cap this at ~8.6 GB (~32% reduction) for ~28% wall-time cost, or `--low-mem-ssd` for a flat ~3.0 GB ceiling at the cost of being a bit slower than `-t 8` overall.
