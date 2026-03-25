#!/usr/bin/env bash
set -euo pipefail

# Benchmark pioSortBed (old & new) vs GNU sort vs bedtools sort.
# Generates random BED files of increasing size, measures wall time and peak RSS,
# and verifies all tools produce output matching LC_ALL=C sort -k1,1 -k2,2n.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OLD="$SCRIPT_DIR/Out/pioSortBed-v1"
NEW="$SCRIPT_DIR/pioSortBed"
TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/pioSortBed-bench.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

SIZES=(10000 100000 1000000 5000000 10000000)
CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr10" "chr11" "chr20" "chrX" "chrY")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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

# Run a command, capture wall time (ms) and peak RSS (kB) via GNU time.
# Usage: run_and_measure <command...>
# For stdin input, set BENCH_STDIN to the file path before calling.
# Sets globals: RESULT_MS, RESULT_KB
TIME_BIN=/usr/bin/time

run_and_measure() {
    local time_out="$TMPDIR/time.tmp"
    if [[ -n "${BENCH_STDIN:-}" ]]; then
        "$TIME_BIN" -f '%e %M' -o "$time_out" "$@" < "$BENCH_STDIN" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp"
    else
        "$TIME_BIN" -f '%e %M' -o "$time_out" "$@" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp"
    fi
    # GNU time outputs: elapsed_seconds peak_rss_kb
    read -r secs kb < "$time_out"
    RESULT_MS=$(awk "BEGIN { printf \"%d\", $secs * 1000 }")
    RESULT_KB=$kb
}

# Like run_and_measure but runs a shell command string (needed for pipelines).
# Usage: run_shell_and_measure "<shell command>"
run_shell_and_measure() {
    local cmd="$1"
    local time_out="$TMPDIR/time.tmp"
    if [[ -n "${BENCH_STDIN:-}" ]]; then
        "$TIME_BIN" -f '%e %M' -o "$time_out" bash -c "$cmd" < "$BENCH_STDIN" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp"
    else
        "$TIME_BIN" -f '%e %M' -o "$time_out" bash -c "$cmd" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp"
    fi
    read -r secs kb < "$time_out"
    RESULT_MS=$(awk "BEGIN { printf \"%d\", $secs * 1000 }")
    RESULT_KB=$kb
}

fmt_reads() {
    local n=$1
    if (( n >= 1000000000 )); then
        awk "BEGIN { printf \"%.0fG\", $n / 1000000000 }"
    elif (( n >= 1000000 )); then
        awk "BEGIN { printf \"%.0fM\", $n / 1000000 }"
    elif (( n >= 1000 )); then
        awk "BEGIN { printf \"%.0fk\", $n / 1000 }"
    else
        echo "$n"
    fi
}

fmt_size() {
    local kb=$1
    if (( kb >= 1048576 )); then
        awk "BEGIN { printf \"%.1f GB\", $kb / 1048576 }"
    elif (( kb >= 1024 )); then
        awk "BEGIN { printf \"%.1f MB\", $kb / 1024 }"
    else
        echo "${kb} kB"
    fi
}

fmt_time() {
    local ms=$1
    if (( ms >= 60000 )); then
        awk "BEGIN { m=int($ms/60000); s=($ms-m*60000)/1000; printf \"%dm%04.1fs\", m, s }"
    elif (( ms >= 10000 )); then
        awk "BEGIN { printf \"%.2fs\", $ms / 1000 }"
    else
        echo "${ms}ms"
    fi
}

# Measure one tool. Sets: RESULT_MS, RESULT_KB
# Usage: bench_one <label> <command...>   — for simple commands
#    or: BENCH_CMD="shell pipeline" bench_one_shell <label>
bench_one() {
    local label="$1"; shift
    cat "$BENCH_FILE" > /dev/null  # warm cache
    run_and_measure "$@"
}

bench_one_shell() {
    local label="$1"; shift
    cat "$BENCH_FILE" > /dev/null
    run_shell_and_measure "$1"
}

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------

if [[ ! -x "$OLD" ]]; then
    echo "Error: old binary not found at $OLD" >&2; exit 1
fi
if [[ ! -x "$NEW" ]]; then
    echo "Error: new binary not found at $NEW (run 'make' first)" >&2; exit 1
fi

HAS_SORT=0
HAS_BEDTOOLS=0
HAS_BEDOPS=0
if command -v sort &>/dev/null; then HAS_SORT=1; fi
if command -v bedtools &>/dev/null; then HAS_BEDTOOLS=1; fi
if command -v sort-bed &>/dev/null; then HAS_BEDOPS=1; fi

echo "pioSortBed benchmark"
echo "===================="
echo "  pioSortBed-v1: $OLD"
echo "  pioSortBed-v2: $NEW"
(( HAS_SORT ))     && echo "  GNU sort:      $(which sort) ($(sort --version | head -1))"
(( HAS_BEDTOOLS )) && echo "  bedtools sort:  $(which bedtools) ($(bedtools --version))"
(( HAS_BEDOPS ))   && echo "  bedops sort-bed: $(which sort-bed) ($(sort-bed --version 2>&1 | head -1))"
echo "  Temp dir:      $TMPDIR"
echo ""

# ---------------------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------------------

# Determine available memory for GNU sort buffer
SORT_BUF="80%"

# CSV output for plotting
CSV_FILE="$SCRIPT_DIR/benchmark_results.csv"
echo "reads,v1_ms,v1_kb,v2_ms,v2_kb,sort_ms,sort_kb,sort1_ms,sort1_kb,bt_ms,bt_kb,bo_ms,bo_kb" > "$CSV_FILE"

SEP="%-10s"
HDR_TIME="  %14s"
HDR_MEM="  %10s"
ROW_TIME="  %14s"
ROW_MEM="  %10s"

# Print header
printf "$SEP" "Reads"
printf "$HDR_TIME$HDR_MEM" "v1 time" "v1 RSS"
printf "$HDR_TIME$HDR_MEM" "v2 time" "v2 RSS"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "sort time" "sort RSS"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "sort1 time" "sort1 RSS"
(( HAS_BEDTOOLS )) && printf "$HDR_TIME$HDR_MEM" "bt time" "bt RSS"
(( HAS_BEDOPS ))   && printf "$HDR_TIME$HDR_MEM" "bo time" "bo RSS"
printf "  %s" "Match"
echo ""

printf "$SEP" "-----"
printf "$HDR_TIME$HDR_MEM" "-------" "------"
printf "$HDR_TIME$HDR_MEM" "-------" "------"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "---------" "--------"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "----------" "---------"
(( HAS_BEDTOOLS )) && printf "$HDR_TIME$HDR_MEM" "-------" "------"
(( HAS_BEDOPS ))   && printf "$HDR_TIME$HDR_MEM" "-------" "------"
printf "  %s" "-----"
echo ""

for n in "${SIZES[@]}"; do
    BENCH_FILE="$TMPDIR/test_${n}.bed"
    BENCH_STDIN=""

    # Generate test data
    generate_bed "$n" "$BENCH_FILE"

    # --- Reference: LC_ALL=C sort -k1,1 -k2,2n (gold standard) ---
    # This defines the expected chromosome order and numeric position order.
    # Not benchmarked in the first pass (timed in the GNU sort slot below).
    LC_ALL=C sort -k1,1 -k2,2n --buffer-size="$SORT_BUF" "$BENCH_FILE" > "$TMPDIR/ref_out.txt"
    # Extract sort keys (chr + beg) from reference for order comparison.
    # Tie-breaking within identical (chr, beg) legitimately varies between tools,
    # so we compare only the key columns for ordering.
    awk -F'\t' '{print $1"\t"$2}' "$TMPDIR/ref_out.txt" > "$TMPDIR/ref_keys.txt"

    # --- pioSortBed v1 (old) ---
    bench_one "v1" "$OLD" "$BENCH_FILE"
    v1_ms=$RESULT_MS; v1_kb=$RESULT_KB
    cp "$TMPDIR/output.tmp" "$TMPDIR/v1_out.txt"

    # --- pioSortBed v2 (new) ---
    bench_one "v2" "$NEW" "$BENCH_FILE"
    v2_ms=$RESULT_MS; v2_kb=$RESULT_KB
    cp "$TMPDIR/output.tmp" "$TMPDIR/v2_out.txt"

    # --- GNU sort multi-threaded (timed) ---
    sort_ms=0; sort_kb=0
    if (( HAS_SORT )); then
        bench_one_shell "sort" "LC_ALL=C sort -k1,1 -k2,2n --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort_ms=$RESULT_MS; sort_kb=$RESULT_KB
    fi

    # --- GNU sort single-threaded (timed) ---
    sort1_ms=0; sort1_kb=0
    if (( HAS_SORT )); then
        bench_one_shell "sort1" "LC_ALL=C sort -k1,1 -k2,2n --parallel=1 --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort1_ms=$RESULT_MS; sort1_kb=$RESULT_KB
    fi

    # --- bedtools sort ---
    bt_ms=0; bt_kb=0
    if (( HAS_BEDTOOLS )); then
        bench_one "bt" bedtools sort -i "$BENCH_FILE"
        bt_ms=$RESULT_MS; bt_kb=$RESULT_KB
        cp "$TMPDIR/output.tmp" "$TMPDIR/bt_out.txt"
    fi

    # --- bedops sort-bed ---
    bo_ms=0; bo_kb=0
    if (( HAS_BEDOPS )); then
        bench_one "bo" sort-bed "$BENCH_FILE"
        bo_ms=$RESULT_MS; bo_kb=$RESULT_KB
        cp "$TMPDIR/output.tmp" "$TMPDIR/bo_out.txt"
    fi

    # --- Verify correctness vs LC_ALL=C sort reference ---
    # 1. Sort-key order: chr+beg columns must match the reference sequence.
    # 2. Line set: fully-sorted output must contain exactly the same lines.
    match="OK"
    for tool in v1 v2 bt bo; do
        outfile="$TMPDIR/${tool}_out.txt"
        [[ -f "$outfile" ]] || continue
        # Check sort-key order matches reference
        awk -F'\t' '{print $1"\t"$2}' "$outfile" > "$TMPDIR/${tool}_keys.txt"
        if ! diff -q "$TMPDIR/ref_keys.txt" "$TMPDIR/${tool}_keys.txt" > /dev/null 2>&1; then
            match="${match}/${tool} ORDER"
            continue
        fi
        # Check same set of lines (fully sort both to normalize tie-breaking)
        LC_ALL=C sort "$outfile" > "$TMPDIR/${tool}_fullsort.txt"
        LC_ALL=C sort "$TMPDIR/ref_out.txt" > "$TMPDIR/ref_fullsort.txt"
        if ! diff -q "$TMPDIR/ref_fullsort.txt" "$TMPDIR/${tool}_fullsort.txt" > /dev/null 2>&1; then
            match="${match}/${tool} LINES"
        fi
    done

    # --- Print row ---
    printf "$SEP" "$(fmt_reads $n)"
    printf "$ROW_TIME$ROW_MEM" "$(fmt_time $v1_ms)" "$(fmt_size $v1_kb)"
    printf "$ROW_TIME$ROW_MEM" "$(fmt_time $v2_ms)" "$(fmt_size $v2_kb)"
    if (( HAS_SORT )); then
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $sort_ms)" "$(fmt_size $sort_kb)"
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $sort1_ms)" "$(fmt_size $sort1_kb)"
    fi
    if (( HAS_BEDTOOLS )); then
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $bt_ms)" "$(fmt_size $bt_kb)"
    fi
    if (( HAS_BEDOPS )); then
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $bo_ms)" "$(fmt_size $bo_kb)"
    fi
    printf "  %s" "$match"
    echo ""

    # Write CSV row
    echo "$n,$v1_ms,$v1_kb,$v2_ms,$v2_kb,$sort_ms,$sort_kb,$sort1_ms,$sort1_kb,$bt_ms,$bt_kb,$bo_ms,$bo_kb" >> "$CSV_FILE"

    # Clean up large files between runs
    rm -f "$BENCH_FILE" "$TMPDIR"/*.txt
done

echo ""

# ---------------------------------------------------------------------------
# Summary: speedup vs v2 (new)
# ---------------------------------------------------------------------------

echo "Speedup relative to pioSortBed-v2 (new):"
echo "(values >1x mean v2 is faster)"
echo ""

# Re-run to collect data for summary (reuse last run if only one size, otherwise regenerate)
printf "%-10s  %8s" "Reads" "vs v1"
(( HAS_SORT ))     && printf "  %8s" "vs sort"
(( HAS_SORT ))     && printf "  %8s" "vs sort1"
(( HAS_BEDTOOLS )) && printf "  %8s" "vs bt"
(( HAS_BEDOPS ))   && printf "  %8s" "vs bo"
echo ""
printf "%-10s  %8s" "-----" "-----"
(( HAS_SORT ))     && printf "  %8s" "-------"
(( HAS_SORT ))     && printf "  %8s" "--------"
(( HAS_BEDTOOLS )) && printf "  %8s" "-----"
(( HAS_BEDOPS ))   && printf "  %8s" "-----"
echo ""

for n in "${SIZES[@]}"; do
    BENCH_FILE="$TMPDIR/test_${n}.bed"
    BENCH_STDIN=""
    generate_bed "$n" "$BENCH_FILE"

    bench_one "v1" "$OLD" "$BENCH_FILE";  v1_ms=$RESULT_MS
    bench_one "v2" "$NEW" "$BENCH_FILE";  v2_ms=$RESULT_MS

    sort_ms=0; sort1_ms=0; bt_ms=0
    if (( HAS_SORT )); then
        bench_one_shell "sort" "LC_ALL=C sort -k1,1 -k2,2n --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort_ms=$RESULT_MS
        bench_one_shell "sort1" "LC_ALL=C sort -k1,1 -k2,2n --parallel=1 --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort1_ms=$RESULT_MS
    fi
    if (( HAS_BEDTOOLS )); then
        bench_one "bt" bedtools sort -i "$BENCH_FILE"
        bt_ms=$RESULT_MS
    fi
    bo_ms=0
    if (( HAS_BEDOPS )); then
        bench_one "bo" sort-bed "$BENCH_FILE"
        bo_ms=$RESULT_MS
    fi

    printf "%-10s" "$(fmt_reads $n)"
    if (( v2_ms > 0 )); then
        printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $v1_ms / $v2_ms }")"
        (( HAS_SORT ))     && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $sort_ms / $v2_ms }")"
        (( HAS_SORT ))     && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $sort1_ms / $v2_ms }")"
        (( HAS_BEDTOOLS )) && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $bt_ms / $v2_ms }")"
        (( HAS_BEDOPS ))   && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $bo_ms / $v2_ms }")"
    else
        printf "  %8s" "inf"
        (( HAS_SORT ))     && printf "  %8s" "inf"
        (( HAS_SORT ))     && printf "  %8s" "inf"
        (( HAS_BEDTOOLS )) && printf "  %8s" "inf"
        (( HAS_BEDOPS ))   && printf "  %8s" "inf"
    fi
    echo ""

    rm -f "$BENCH_FILE"
done

echo ""

# ---------------------------------------------------------------------------
# Stdin vs file comparison (new version only)
# ---------------------------------------------------------------------------

echo "Stdin vs file I/O (pioSortBed-v2, 1M reads)"
echo "----------------------------------------------"
BENCH_FILE="$TMPDIR/test_stdin.bed"
generate_bed 1000000 "$BENCH_FILE"

BENCH_STDIN=""
bench_one "file" "$NEW" "$BENCH_FILE"
file_ms=$RESULT_MS; file_kb=$RESULT_KB

BENCH_STDIN="$BENCH_FILE"
bench_one "stdin" "$NEW" "-"
stdin_ms=$RESULT_MS; stdin_kb=$RESULT_KB
BENCH_STDIN=""

printf "  File:  %s, %s peak RSS\n" "$(fmt_time $file_ms)" "$(fmt_size $file_kb)"
printf "  Stdin: %s, %s peak RSS\n" "$(fmt_time $stdin_ms)" "$(fmt_size $stdin_kb)"

echo ""

# ---------------------------------------------------------------------------
# Generate plot with gnuplot
# ---------------------------------------------------------------------------

PLOT_FILE="$SCRIPT_DIR/benchmark_plot.png"

if command -v gnuplot &>/dev/null; then
    echo "Generating plot: $PLOT_FILE"
    GNUPLOT_SCRIPT="$TMPDIR/plot.gp"
    cat > "$GNUPLOT_SCRIPT" <<'GPEOF'
set terminal pngcairo size 900,500 enhanced font 'Arial,11'
set output 'PLOT_PLACEHOLDER'

set datafile separator ','

set title "pioSortBed Benchmark" font ',14'
set xlabel 'Number of reads'
set ylabel 'Wall time (seconds)' textcolor rgb '#333333'
set y2label 'Peak RSS (MB)' textcolor rgb '#999999'

set logscale x 10
set logscale y 10
set logscale y2 10
set format x '%.0s%c'
set format y '%.2g'
set format y2 '%.0f'

set xtics nomirror
set ytics nomirror textcolor rgb '#333333'
set y2tics nomirror textcolor rgb '#999999'

set key top left font ',10' spacing 1.2
set grid xtics ytics lt 0 lw 0.5 lc rgb '#dddddd'
set border 11

set style line 1 lc rgb '#e41a1c' lw 2.2 pt 7  ps 1.0
set style line 2 lc rgb '#377eb8' lw 2.2 pt 7  ps 1.0
set style line 3 lc rgb '#4daf4a' lw 2.2 pt 7  ps 1.0
set style line 4 lc rgb '#ff7f00' lw 2.2 pt 7  ps 1.0
set style line 5 lc rgb '#984ea3' lw 2.2 pt 7  ps 1.0
set style line 12 lc rgb '#377eb8' lw 1.4 pt 5 ps 0.8 dt 4
set style line 13 lc rgb '#4daf4a' lw 1.4 pt 5 ps 0.8 dt 4
set style line 14 lc rgb '#ff7f00' lw 1.4 pt 5 ps 0.8 dt 4
set style line 15 lc rgb '#984ea3' lw 1.4 pt 5 ps 0.8 dt 4

plot 'CSV_PLACEHOLDER' \
    skip 1 using 1:($2/1000.0) axes x1y1 with linespoints ls 1 title 'pioSortBed v1', \
 '' skip 1 using 1:($4/1000.0) axes x1y1 with linespoints ls 2 title 'pioSortBed v2', \
 '' skip 1 using 1:($6/1000.0) axes x1y1 with linespoints ls 3 title 'GNU sort', \
 '' skip 1 using 1:($10/1000.0) axes x1y1 with linespoints ls 4 title 'bedtools sort', \
 '' skip 1 using 1:($12/1000.0) axes x1y1 with linespoints ls 5 title 'bedops sort-bed', \
 '' skip 1 using 1:($5/1024.0)  axes x1y2 with linespoints ls 12 title 'v2 mem', \
 '' skip 1 using 1:($7/1024.0)  axes x1y2 with linespoints ls 13 title 'sort mem', \
 '' skip 1 using 1:($11/1024.0) axes x1y2 with linespoints ls 14 title 'bedtools mem', \
 '' skip 1 using 1:($13/1024.0) axes x1y2 with linespoints ls 15 title 'bedops mem'
GPEOF
    sed -i "s|PLOT_PLACEHOLDER|$PLOT_FILE|g; s|CSV_PLACEHOLDER|$CSV_FILE|g" "$GNUPLOT_SCRIPT"
    gnuplot "$GNUPLOT_SCRIPT" && echo "Plot saved to $PLOT_FILE" || echo "gnuplot failed"
else
    echo "gnuplot not found — skipping plot generation"
fi

echo ""
echo "Done."
