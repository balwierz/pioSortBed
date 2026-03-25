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

SIZES=(10000 100000 1000000 5000000 10000000 50000000)
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

# Run a command, capture wall time (ms) and peak RSS (kB) via /proc.
# Usage: run_and_measure <command...>
# For stdin input, set BENCH_STDIN to the file path before calling.
# Sets globals: PEAK_KB
run_and_measure() {
    if [[ -n "${BENCH_STDIN:-}" ]]; then
        "$@" < "$BENCH_STDIN" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp" &
    else
        "$@" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp" &
    fi
    local pid=$!

    local peak=0
    while kill -0 "$pid" 2>/dev/null; do
        local rss
        rss=$(awk '/^VmHWM:/ {print $2}' "/proc/$pid/status" 2>/dev/null) || true
        if [[ -n "$rss" && "$rss" -gt "$peak" ]]; then
            peak=$rss
        fi
        local rss_cur
        rss_cur=$(awk '/^VmRSS:/ {print $2}' "/proc/$pid/status" 2>/dev/null) || true
        if [[ -n "$rss_cur" && "$rss_cur" -gt "$peak" ]]; then
            peak=$rss_cur
        fi
        sleep 0.01
    done
    wait "$pid" || true

    PEAK_KB=$peak
}

# Like run_and_measure but runs a shell command string (needed for pipelines).
# Usage: run_shell_and_measure "<shell command>"
run_shell_and_measure() {
    local cmd="$1"
    if [[ -n "${BENCH_STDIN:-}" ]]; then
        bash -c "$cmd" < "$BENCH_STDIN" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp" &
    else
        bash -c "$cmd" > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp" &
    fi
    local pid=$!

    local peak=0
    while kill -0 "$pid" 2>/dev/null; do
        local rss
        rss=$(awk '/^VmHWM:/ {print $2}' "/proc/$pid/status" 2>/dev/null) || true
        if [[ -n "$rss" && "$rss" -gt "$peak" ]]; then
            peak=$rss
        fi
        local rss_cur
        rss_cur=$(awk '/^VmRSS:/ {print $2}' "/proc/$pid/status" 2>/dev/null) || true
        if [[ -n "$rss_cur" && "$rss_cur" -gt "$peak" ]]; then
            peak=$rss_cur
        fi
        sleep 0.01
    done
    wait "$pid" || true

    PEAK_KB=$peak
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
    local t0=$(($(date +%s%N)))
    run_and_measure "$@"
    local t1=$(($(date +%s%N)))
    RESULT_MS=$(( (t1 - t0) / 1000000 ))
    RESULT_KB=$PEAK_KB
}

bench_one_shell() {
    local label="$1"; shift
    cat "$BENCH_FILE" > /dev/null
    local t0=$(($(date +%s%N)))
    run_shell_and_measure "$1"
    local t1=$(($(date +%s%N)))
    RESULT_MS=$(( (t1 - t0) / 1000000 ))
    RESULT_KB=$PEAK_KB
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
if command -v sort &>/dev/null; then HAS_SORT=1; fi
if command -v bedtools &>/dev/null; then HAS_BEDTOOLS=1; fi

echo "pioSortBed benchmark"
echo "===================="
echo "  pioSortBed-v1: $OLD"
echo "  pioSortBed-v2: $NEW"
(( HAS_SORT ))     && echo "  GNU sort:      $(which sort) ($(sort --version | head -1))"
(( HAS_BEDTOOLS )) && echo "  bedtools sort:  $(which bedtools) ($(bedtools --version))"
echo "  Temp dir:      $TMPDIR"
echo ""

# ---------------------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------------------

# Determine available memory for GNU sort buffer
SORT_BUF="80%"

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
printf "  %s" "Match"
echo ""

printf "$SEP" "-----"
printf "$HDR_TIME$HDR_MEM" "-------" "------"
printf "$HDR_TIME$HDR_MEM" "-------" "------"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "---------" "--------"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "----------" "---------"
(( HAS_BEDTOOLS )) && printf "$HDR_TIME$HDR_MEM" "-------" "------"
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

    # --- Verify correctness vs LC_ALL=C sort reference ---
    # 1. Sort-key order: chr+beg columns must match the reference sequence.
    # 2. Line set: fully-sorted output must contain exactly the same lines.
    match="OK"
    for tool in v1 v2 bt; do
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
    printf "$SEP" "$n"
    printf "$ROW_TIME$ROW_MEM" "$(fmt_time $v1_ms)" "$(fmt_size $v1_kb)"
    printf "$ROW_TIME$ROW_MEM" "$(fmt_time $v2_ms)" "$(fmt_size $v2_kb)"
    if (( HAS_SORT )); then
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $sort_ms)" "$(fmt_size $sort_kb)"
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $sort1_ms)" "$(fmt_size $sort1_kb)"
    fi
    if (( HAS_BEDTOOLS )); then
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $bt_ms)" "$(fmt_size $bt_kb)"
    fi
    printf "  %s" "$match"
    echo ""

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
echo ""
printf "%-10s  %8s" "-----" "-----"
(( HAS_SORT ))     && printf "  %8s" "-------"
(( HAS_SORT ))     && printf "  %8s" "--------"
(( HAS_BEDTOOLS )) && printf "  %8s" "-----"
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

    printf "%-10s" "$n"
    if (( v2_ms > 0 )); then
        printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $v1_ms / $v2_ms }")"
        (( HAS_SORT ))     && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $sort_ms / $v2_ms }")"
        (( HAS_SORT ))     && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $sort1_ms / $v2_ms }")"
        (( HAS_BEDTOOLS )) && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $bt_ms / $v2_ms }")"
    else
        printf "  %8s" "inf"
        (( HAS_SORT ))     && printf "  %8s" "inf"
        (( HAS_SORT ))     && printf "  %8s" "inf"
        (( HAS_BEDTOOLS )) && printf "  %8s" "inf"
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
echo "Done."
