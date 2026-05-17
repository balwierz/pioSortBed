#!/usr/bin/env bash
set -euo pipefail

# Benchmark Loci (Python; polars + pandas sort backends) across the same
# input sizes as benchmark.sh, for direct apples-to-apples scaling
# comparison against pioSortBed / GNU sort / bedtools / bedops.
#
# Why a dedicated script:
# - Loci adds ~0.3 s of Python interpreter + import startup per invocation,
#   which dominates at the small end of the size sweep and would distort
#   plot_readme.gp's log-y axis. Keeping Loci in its own CSV/plot lets
#   readers see the scaling curve without rescaling the main plot.
# - Avoids modifying benchmark.sh's CSV schema (consumed by plot_readme.gp,
#   archived under benchmark/history/, and post-processed by an awk script
#   that hard-codes column positions). Adding columns there is a chain of
#   coupled changes; this script just emits its own CSV.
#
# Sizes match benchmark.sh by default; override with SIZES="10000 100000".
# Loci wall time is reported as the wall of the whole script invocation,
# including Python import + sort + write. This is the honest user-facing
# cost.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/loci-bench.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

LOCI_PY="${LOCI_PY:-/home/piotr/Environments/Python3.14/bin/python}"
LOCI_SORT_SCRIPT="${LOCI_SORT_SCRIPT:-/home/piotr/Sources/Loci1/examples/sort_bed.py}"
TIME_BIN=/usr/bin/time

if [[ ! -x "$LOCI_PY" ]]; then
    echo "Error: Loci python not found at $LOCI_PY" >&2
    echo "Override with LOCI_PY=/path/to/python" >&2
    exit 1
fi
if [[ ! -f "$LOCI_SORT_SCRIPT" ]]; then
    echo "Error: sort_bed.py not found at $LOCI_SORT_SCRIPT" >&2
    echo "Override with LOCI_SORT_SCRIPT=/path/to/sort_bed.py" >&2
    exit 1
fi
LOCI_VER=$("$LOCI_PY" -c 'import loci; print(loci.__version__)' 2>/dev/null \
           || { echo "Error: 'import loci' failed in $LOCI_PY" >&2; exit 1; })

# Size sweep — same as benchmark.sh by default. Override via env, e.g.
#   SIZES="10000 100000 1000000" bash bench_loci.sh
if [[ -n "${SIZES:-}" ]]; then
    read -r -a SIZES <<< "$SIZES"
else
    SIZES=(10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000 200000000)
fi

# Above this size, drop -t 1 polars/pandas runs (they take many minutes at
# 100M+ and pandas' lex-sort goes super-linear). Override via HUGE_SKIP.
HUGE_SKIP="${HUGE_SKIP:-100000000}"

CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr10" "chr11" "chr20" "chrX" "chrY")

generate_bed() {
    local n=$1 out=$2
    awk -v n="$n" -v nchr=${#CHROMS[@]} '
    BEGIN {
        split("'"$(IFS=,; echo "${CHROMS[*]}")"'", chrs, ",")
        srand(42)
        for (i = 0; i < n; i++) {
            c = chrs[int(rand() * nchr) + 1]
            beg = int(rand() * 249000000)
            end = beg + int(rand() * 1000) + 1
            printf "%s\t%d\t%d\tread%d\t%d\t%s\n", c, beg, end, i, int(rand()*1000), (rand()>0.5 ? "+" : "-")
        }
    }' > "$out"
}

run_loci() {
    # run_loci <backend> <cores> <input> -> sets RESULT_MS, RESULT_KB
    local backend="$1" cores="$2" input="$3"
    local time_out="$TMPDIR/time.tmp"
    "$TIME_BIN" -f '%e %M' -o "$time_out" \
        "$LOCI_PY" "$LOCI_SORT_SCRIPT" "$input" -c "$cores" -b "$backend" \
        > /dev/null 2>"$TMPDIR/stderr.tmp" || {
        RESULT_MS="NA"; RESULT_KB="NA"
        return
    }
    read -r secs kb < "$time_out"
    RESULT_MS=$(awk "BEGIN { printf \"%d\", $secs * 1000 }")
    RESULT_KB=$kb
}

fmt_reads() {
    local n=$1
    if (( n >= 1000000 )); then awk "BEGIN { printf \"%.0fM\", $n / 1000000 }"
    elif (( n >= 1000 )); then awk "BEGIN { printf \"%.0fk\", $n / 1000 }"
    else echo "$n"
    fi
}
fmt_ms() { local ms=$1; [[ "$ms" = NA ]] && { echo "—"; return; }; awk "BEGIN { printf \"%.2fs\", $ms/1000.0 }"; }
fmt_mb() { local kb=$1; [[ "$kb" = NA ]] && { echo "—"; return; }; awk "BEGIN { printf \"%.0fMB\", $kb/1024.0 }"; }

echo "==============================================================="
echo "Loci sort-time scaling benchmark"
echo "  loci   : $LOCI_PY (loci $LOCI_VER)"
echo "  script : $LOCI_SORT_SCRIPT"
echo "  sizes  : ${SIZES[*]}"
echo "  HUGE_SKIP (drop -t 1 above) : $HUGE_SKIP"
echo "==============================================================="

CSV="$SCRIPT_DIR/bench_loci.csv"
echo "reads,polars1_ms,polars1_kb,polars8_ms,polars8_kb,pandas1_ms,pandas1_kb,pandas8_ms,pandas8_kb" > "$CSV"

printf "\n%-8s  %12s %10s  %12s %10s  %12s %10s  %12s %10s\n" \
    "Reads" "polars-1t" "RSS" "polars-8t" "RSS" "pandas-1t" "RSS" "pandas-8t" "RSS"
printf '%s\n' "-----------------------------------------------------------------------------------------------------------"

for n in "${SIZES[@]}"; do
    BENCH_FILE="$TMPDIR/test_${n}.bed"
    generate_bed "$n" "$BENCH_FILE"

    if (( n >= HUGE_SKIP )); then
        polars1_ms=NA; polars1_kb=NA
        pandas1_ms=NA; pandas1_kb=NA
    else
        run_loci polars 1 "$BENCH_FILE"; polars1_ms=$RESULT_MS; polars1_kb=$RESULT_KB
        run_loci pandas 1 "$BENCH_FILE"; pandas1_ms=$RESULT_MS; pandas1_kb=$RESULT_KB
    fi
    run_loci polars 8 "$BENCH_FILE"; polars8_ms=$RESULT_MS; polars8_kb=$RESULT_KB
    run_loci pandas 8 "$BENCH_FILE"; pandas8_ms=$RESULT_MS; pandas8_kb=$RESULT_KB

    printf "%-8s  %12s %10s  %12s %10s  %12s %10s  %12s %10s\n" \
        "$(fmt_reads $n)" \
        "$(fmt_ms $polars1_ms)" "$(fmt_mb $polars1_kb)" \
        "$(fmt_ms $polars8_ms)" "$(fmt_mb $polars8_kb)" \
        "$(fmt_ms $pandas1_ms)" "$(fmt_mb $pandas1_kb)" \
        "$(fmt_ms $pandas8_ms)" "$(fmt_mb $pandas8_kb)"

    echo "$n,$polars1_ms,$polars1_kb,$polars8_ms,$polars8_kb,$pandas1_ms,$pandas1_kb,$pandas8_ms,$pandas8_kb" >> "$CSV"

    rm -f "$BENCH_FILE"
done

echo ""
echo "Wrote $CSV"
echo "Done."
