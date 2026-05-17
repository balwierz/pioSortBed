#!/usr/bin/env bash
# Query-performance benchmark: region queries on three indexed forms of
# the same sorted BED — bgzip+tabix, LociSSD (no index), LociSSD + interval
# index — using four readers — tabix CLI / polars / DuckDB / Loci.
#
# Inputs:
#   $1: workspace dir holding regions.bed.gz / regions.bed.gz.tbi /
#       regions.lociss / regions_idx.lociss. Default /var/tmp/qbench.
#   $2: chromosome lengths file (TSV: chrom<TAB>length). Default
#       benchmark/grch38_chromsizes.tsv (generated below if missing).
#
# Outputs:
#   benchmark/bench_query.csv     — one row per (form × reader × size × cache)
#   benchmark/bench_query_queries.tsv — the random query workload, for re-use
#
set -euo pipefail

cd "$(dirname "$0")/.."

QDIR="${1:-/var/tmp/qbench}"
CHROMSIZES="${2:-benchmark/grch38_chromsizes.tsv}"
PY=/home/piotr/Environments/Python3.14/bin/python
DRIVER="$(pwd)/benchmark/bench_query_drivers.py"
CSV="$(pwd)/benchmark/bench_query.csv"
QFILE="$(pwd)/benchmark/bench_query_queries.tsv"

for f in "$QDIR/regions.bed.gz" "$QDIR/regions.bed.gz.tbi" \
         "$QDIR/regions.lociss" "$QDIR/regions_idx.lociss" \
         "$QDIR/regions10m.bed.gz" "$QDIR/regions10m.bed.gz.tbi" \
         "$QDIR/regions10m.lociss" "$QDIR/regions10m_idx.lociss"; do
    [[ -f "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }
done

# ----------------------------------------------------------------------
# Generate chromosome-lengths table from the LociSSD manifest if missing.
# (Each query needs a chromosome with a known max length so we can pick a
# valid random start.)
# ----------------------------------------------------------------------
if [[ ! -s "$CHROMSIZES" ]]; then
    echo "Generating $CHROMSIZES from LociSSD manifest..."
    "$PY" - <<PY
import json, pyarrow.parquet as pq, sys
md = pq.read_metadata("$QDIR/regions.lociss")
manifest = json.loads(md.metadata[b"lociSSD_manifest"].decode("utf-8"))
with open("$CHROMSIZES", "w") as fh:
    for c in manifest["chromosomes"]:
        # max_end is the largest End across all rows for this chromosome —
        # use it as the chromosome's effective span for query generation.
        fh.write(f"{c['name']}\t{c['max_end']}\n")
PY
fi

# ----------------------------------------------------------------------
# Generate the 1000-query workload (333 × 3 region sizes) — once, cached.
# ----------------------------------------------------------------------
if [[ ! -s "$QFILE" ]]; then
    echo "Generating query workload at $QFILE..."
    "$PY" - <<PY
import random
random.seed(42)
chrs = []
for line in open("$CHROMSIZES"):
    name, length = line.rstrip().split("\t")
    chrs.append((name, int(length)))

# 333 of each: 1 kbp, 100 kbp, 10 Mbp
sizes = [1_000] * 333 + [100_000] * 333 + [10_000_000] * 334
random.shuffle(sizes)

with open("$QFILE", "w") as fh:
    for size in sizes:
        # pick a chrom large enough to host this region
        candidates = [(n, L) for (n, L) in chrs if L > size]
        if not candidates:
            continue
        name, L = random.choice(candidates)
        start = random.randint(0, L - size - 1)
        end = start + size
        fh.write(f"{name}\t{start}\t{end}\t{size}\n")
PY
fi
QCOUNT=$(wc -l < "$QFILE")
echo "Workload: $QCOUNT queries"

# ----------------------------------------------------------------------
# Helper: drop FS cache for a path (best-effort, no root).
# ----------------------------------------------------------------------
drop_fcache() {
    sync
    for f in "$@"; do
        [[ -f "$f" ]] && dd if="$f" of=/dev/null bs=1M iflag=nocache 2>/dev/null \
            || true
    done
}

# ----------------------------------------------------------------------
# Run one (reader × form × size × cache) cell. Splits the workload by
# region size and feeds only matching queries to the driver. Records:
#   median_us, p99_us, qcount.
# ----------------------------------------------------------------------
run_cell() {
    local reader="$1" path="$2" size="$3" cache_state="$4" extra_arg="${5:-}"
    local lat="$QDIR/lat_${reader}_${size}_${cache_state}.txt"
    rm -f "$lat"

    # Filter queries to this size bucket.
    local size_queries
    size_queries=$(mktemp)
    awk -F'\t' -v sz="$size" '$4 == sz {print $1"\t"$2"\t"$3}' "$QFILE" \
        > "$size_queries"
    local n
    n=$(wc -l < "$size_queries")

    if [[ "$cache_state" == "cold" ]]; then
        # Cold-cache: drop FS cache before the run. We don't drop between
        # queries (would need root to be reliable); the cold metric here
        # captures "open + read first time".
        drop_fcache "$path" "${path}.tbi"
    fi

    "$PY" "$DRIVER" "$reader" "$path" $extra_arg < "$size_queries" > "$lat"

    # Compute median + p99 + total time (for throughput).
    local stats
    stats=$("$PY" - <<PY
import statistics, sys
xs = [int(l) for l in open("$lat") if l.strip()]
xs.sort()
n = len(xs)
if n == 0:
    print("0\t0\t0\t0")
else:
    median = xs[n//2]
    p99    = xs[min(int(0.99*n), n-1)]
    total  = sum(xs)
    print(f"{median}\t{p99}\t{total}\t{n}")
PY
)
    rm -f "$size_queries"
    local median_us p99_us total_us cells
    IFS=$'\t' read -r median_us p99_us total_us cells <<< "$stats"
    printf "  %-8s %-8s %-12s cache=%-4s n=%-4s median=%-8s p99=%-8s tot=%s us\n" \
        "$reader" "$size" "$(basename "$path")" "$cache_state" "$cells" \
        "$median_us" "$p99_us" "$total_us"
    # CSV row
    local throughput
    throughput=$(awk -v t="$total_us" -v n="$cells" \
        'BEGIN{if(t>0)printf "%.1f", n*1e6/t; else print "0"}')
    echo "$reader,$(basename "$path"),$size,$cache_state,$cells,$median_us,$p99_us,$throughput" >> "$CSV"
}

# ----------------------------------------------------------------------
# Run the sweep.
# ----------------------------------------------------------------------
echo "reader,form,region_size,cache_state,n_queries,median_us,p99_us,qps" > "$CSV"

# NOTE: at 100M records the embedded interval index is multi-GB and
# exceeds pyarrow's default thrift footer size limit, so loci can't
# open the indexed file. This is documented in FORMAT_SPEC.md §6.5.1
# (the index lives in Parquet footer KV metadata; large index footers
# are hostile to pyarrow's defaults). The 100M cell skips the
# "loci+index" row; we report the with-index advantage at 10M instead,
# where the index fits comfortably (226 MB total file, ~67 MB of
# index in the footer).
for cache in warm cold; do
    for size in 1000 100000 10000000; do
        echo ""
        echo "=== cache=$cache size=${size} | 100M-row workload ==="
        run_cell tabix  "$QDIR/regions.bed.gz"      "$size" "$cache"
        run_cell polars "$QDIR/regions.lociss"      "$size" "$cache"
        run_cell duckdb "$QDIR/regions.lociss"      "$size" "$cache"
        run_cell loci   "$QDIR/regions.lociss"      "$size" "$cache" ""

        echo ""
        echo "=== cache=$cache size=${size} | 10M-row workload (with --lociss-index) ==="
        run_cell tabix  "$QDIR/regions10m.bed.gz"     "$size" "$cache"
        run_cell loci   "$QDIR/regions10m.lociss"     "$size" "$cache" ""
        run_cell loci   "$QDIR/regions10m_idx.lociss" "$size" "$cache" --use-index
    done
done

echo ""
echo "Wrote: $CSV"
