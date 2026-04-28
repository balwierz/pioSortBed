#!/usr/bin/env bash
set -euo pipefail

# Benchmark pioSortBed vs GNU sort vs bedtools sort vs bedops sort-bed.
# Generates random BED files of increasing size, measures wall time and peak RSS,
# and verifies all tools produce output matching LC_ALL=C sort -k1,1 -k2,2n.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIO="$REPO_DIR/pioSortBed"
TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/pioSortBed-bench.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

# Parse options
VERIFY=0
for arg in "$@"; do
    case "$arg" in
        --verify)    VERIFY=1 ;;
        --no-verify) VERIFY=0 ;;
    esac
done

SIZES=(100000 1000000 5000000 10000000 50000000)
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
        awk "BEGIN { m=int($ms/60000); s=($ms-m*60000)/1000; printf \"%dmin%04.1fs\", m, s }"
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

if [[ ! -x "$PIO" ]]; then
    echo "Error: pioSortBed binary not found at $PIO (run 'make' first)" >&2; exit 1
fi

HAS_SORT=0
HAS_BEDTOOLS=0
HAS_BEDOPS=0
if command -v sort &>/dev/null; then HAS_SORT=1; fi
if command -v bedtools &>/dev/null; then HAS_BEDTOOLS=1; fi
if command -v sort-bed &>/dev/null; then HAS_BEDOPS=1; fi

echo "pioSortBed benchmark"
echo "===================="
echo "  pioSortBed:    $PIO"
(( HAS_SORT ))     && echo "  GNU sort:      $(which sort) ($(sort --version | head -1))"
(( HAS_BEDTOOLS )) && echo "  bedtools sort:  $(which bedtools) ($(bedtools --version))"
if (( HAS_BEDOPS )); then
    bedops_ver=$(sort-bed --version 2>&1 | awk -F'version:' '/version:/ {gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}')
    echo "  bedops sort-bed: $(which sort-bed) (version: $bedops_ver)"
fi
echo "  Temp dir:      $TMPDIR"
echo ""

# ---------------------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------------------

# Determine available memory for GNU sort buffer
SORT_BUF="80%"

# CSV output for plotting
CSV_FILE="$SCRIPT_DIR/benchmark_results.csv"
echo "reads,pio1_ms,pio1_kb,pio8_ms,pio8_kb,pio_lm_ms,pio_lm_kb,sort1_ms,sort1_kb,sort8_ms,sort8_kb,bt_ms,bt_kb,bo_ms,bo_kb" > "$CSV_FILE"

SEP="%-10s"
HDR_TIME="  %14s"
HDR_MEM="  %10s"
ROW_TIME="  %14s"
ROW_MEM="  %10s"

# Print header
printf "$SEP" "Reads"
printf "$HDR_TIME$HDR_MEM" "pio-1t" "RSS"
printf "$HDR_TIME$HDR_MEM" "pio-8t" "RSS"
printf "$HDR_TIME$HDR_MEM" "pio-lm" "RSS"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "sort-1t" "RSS"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "sort-8t" "RSS"
(( HAS_BEDTOOLS )) && printf "$HDR_TIME$HDR_MEM" "bedtools" "RSS"
(( HAS_BEDOPS ))   && printf "$HDR_TIME$HDR_MEM" "bedops" "RSS"
printf "  %s" "Match"
echo ""

printf "$SEP" "-----"
printf "$HDR_TIME$HDR_MEM" "------" "---"
printf "$HDR_TIME$HDR_MEM" "------" "---"
printf "$HDR_TIME$HDR_MEM" "------" "---"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "-------" "---"
(( HAS_SORT ))     && printf "$HDR_TIME$HDR_MEM" "-------" "---"
(( HAS_BEDTOOLS )) && printf "$HDR_TIME$HDR_MEM" "--------" "---"
(( HAS_BEDOPS ))   && printf "$HDR_TIME$HDR_MEM" "------" "---"
printf "  %s" "-----"
echo ""

for n in "${SIZES[@]}"; do
    BENCH_FILE="$TMPDIR/test_${n}.bed"
    BENCH_STDIN=""

    # Generate test data
    generate_bed "$n" "$BENCH_FILE"

    # --- Reference: LC_ALL=C sort -k1,1 -k2,2n (gold standard) ---
    if (( VERIFY )); then
        LC_ALL=C sort -k1,1 -k2,2n --buffer-size="$SORT_BUF" "$BENCH_FILE" > "$TMPDIR/ref_out.txt"
        awk -F'\t' '{print $1"\t"$2}' "$TMPDIR/ref_out.txt" > "$TMPDIR/ref_keys.txt"
    fi


    # --- pioSortBed single-threaded ---
    bench_one "pio1" "$PIO" -t 1 "$BENCH_FILE"
    pio1_ms=$RESULT_MS; pio1_kb=$RESULT_KB
    cp "$TMPDIR/output.tmp" "$TMPDIR/pio1_out.txt"

    # --- pioSortBed 8-threaded (parallelizes at all sizes since v2.1.0) ---
    bench_one "pio8" "$PIO" -t 8 "$BENCH_FILE"
    pio8_ms=$RESULT_MS; pio8_kb=$RESULT_KB
    cp "$TMPDIR/output.tmp" "$TMPDIR/pio8_out.txt"

    # --- pioSortBed low-memory SSD mode ---
    bench_one "pio-lm" "$PIO" --low-mem-ssd "$BENCH_FILE"
    pio_lm_ms=$RESULT_MS; pio_lm_kb=$RESULT_KB
    cp "$TMPDIR/output.tmp" "$TMPDIR/pio_lm_out.txt"

    # --- GNU sort single-threaded ---
    sort1_ms=0; sort1_kb=0
    if (( HAS_SORT )); then
        bench_one_shell "sort1" "LC_ALL=C sort -k1,1 -k2,2n --parallel=1 --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort1_ms=$RESULT_MS; sort1_kb=$RESULT_KB
    fi

    # --- GNU sort 8-threaded ---
    sort8_ms=0; sort8_kb=0
    if (( HAS_SORT )); then
        bench_one_shell "sort8" "LC_ALL=C sort -k1,1 -k2,2n --parallel=8 --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort8_ms=$RESULT_MS; sort8_kb=$RESULT_KB
    fi

    # --- bedtools sort ---
    bt_ms=0; bt_kb=0
    if (( HAS_BEDTOOLS )) && (( n < 100000000 )); then
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
    match="-"
    if (( VERIFY )); then
        match="OK"
        for tool in pio1 pio8 pio_lm bt bo; do
            outfile="$TMPDIR/${tool}_out.txt"
            [[ -f "$outfile" ]] || continue
            awk -F'\t' '{print $1"\t"$2}' "$outfile" > "$TMPDIR/${tool}_keys.txt"
            if ! diff -q "$TMPDIR/ref_keys.txt" "$TMPDIR/${tool}_keys.txt" > /dev/null 2>&1; then
                match="${match}/${tool} ORDER"
                continue
            fi
            LC_ALL=C sort "$outfile" > "$TMPDIR/${tool}_fullsort.txt"
            LC_ALL=C sort "$TMPDIR/ref_out.txt" > "$TMPDIR/ref_fullsort.txt"
            if ! diff -q "$TMPDIR/ref_fullsort.txt" "$TMPDIR/${tool}_fullsort.txt" > /dev/null 2>&1; then
                match="${match}/${tool} LINES"
            fi
        done
    fi

    # --- Print row ---
    printf "$SEP" "$(fmt_reads $n)"
    printf "$ROW_TIME$ROW_MEM" "$(fmt_time $pio1_ms)" "$(fmt_size $pio1_kb)"
    printf "$ROW_TIME$ROW_MEM" "$(fmt_time $pio8_ms)" "$(fmt_size $pio8_kb)"
    printf "$ROW_TIME$ROW_MEM" "$(fmt_time $pio_lm_ms)" "$(fmt_size $pio_lm_kb)"
    if (( HAS_SORT )); then
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $sort1_ms)" "$(fmt_size $sort1_kb)"
        printf "$ROW_TIME$ROW_MEM" "$(fmt_time $sort8_ms)" "$(fmt_size $sort8_kb)"
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
    echo "$n,$pio1_ms,$pio1_kb,$pio8_ms,$pio8_kb,$pio_lm_ms,$pio_lm_kb,$sort1_ms,$sort1_kb,$sort8_ms,$sort8_kb,$bt_ms,$bt_kb,$bo_ms,$bo_kb" >> "$CSV_FILE"

    # Clean up large files between runs
    rm -f "$BENCH_FILE" "$TMPDIR"/*.txt
done

echo ""

# ---------------------------------------------------------------------------
# Summary: speedup vs v2 (new)
# ---------------------------------------------------------------------------

echo "Speedup relative to pioSortBed 8-thread:"
echo "(values >1x mean pioSortBed-8t is faster)"
echo ""

printf "%-10s  %8s" "Reads" "vs pio1"
printf "  %8s" "vs pio-lm"
(( HAS_SORT ))     && printf "  %8s" "vs sort1"
(( HAS_SORT ))     && printf "  %8s" "vs sort8"
(( HAS_BEDTOOLS )) && printf "  %8s" "vs bt"
(( HAS_BEDOPS ))   && printf "  %8s" "vs bo"
echo ""
printf "%-10s  %8s" "-----" "-------"
printf "  %8s" "---------"
(( HAS_SORT ))     && printf "  %8s" "-------"
(( HAS_SORT ))     && printf "  %8s" "-------"
(( HAS_BEDTOOLS )) && printf "  %8s" "-----"
(( HAS_BEDOPS ))   && printf "  %8s" "-----"
echo ""

for n in "${SIZES[@]}"; do
    BENCH_FILE="$TMPDIR/test_${n}.bed"
    BENCH_STDIN=""
    generate_bed "$n" "$BENCH_FILE"

    bench_one "pio1" "$PIO" -t 1 "$BENCH_FILE";  pio1_ms=$RESULT_MS
    bench_one "pio8" "$PIO" -t 8 "$BENCH_FILE";  pio8_ms=$RESULT_MS
    bench_one "pio-lm" "$PIO" --low-mem-ssd "$BENCH_FILE";  pio_lm_ms=$RESULT_MS

    sort1_ms=0; sort8_ms=0; bt_ms=0; bo_ms=0
    if (( HAS_SORT )); then
        bench_one_shell "sort1" "LC_ALL=C sort -k1,1 -k2,2n --parallel=1 --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort1_ms=$RESULT_MS
        bench_one_shell "sort8" "LC_ALL=C sort -k1,1 -k2,2n --parallel=8 --buffer-size=$SORT_BUF '$BENCH_FILE'"
        sort8_ms=$RESULT_MS
    fi
    if (( HAS_BEDTOOLS )) && (( n < 100000000 )); then
        bench_one "bt" bedtools sort -i "$BENCH_FILE"
        bt_ms=$RESULT_MS
    fi
    if (( HAS_BEDOPS )); then
        bench_one "bo" sort-bed "$BENCH_FILE"
        bo_ms=$RESULT_MS
    fi

    printf "%-10s" "$(fmt_reads $n)"
    if (( pio8_ms > 0 )); then
        printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $pio1_ms / $pio8_ms }")"
        printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $pio_lm_ms / $pio8_ms }")"
        (( HAS_SORT ))     && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $sort1_ms / $pio8_ms }")"
        (( HAS_SORT ))     && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $sort8_ms / $pio8_ms }")"
        (( HAS_BEDTOOLS )) && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $bt_ms / $pio8_ms }")"
        (( HAS_BEDOPS ))   && printf "  %8s" "$(awk "BEGIN { printf \"%.2fx\", $bo_ms / $pio8_ms }")"
    else
        printf "  %8s" "inf"
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

echo "Stdin vs file I/O (pioSortBed, 1M reads)"
echo "-------------------------------------------"
BENCH_FILE="$TMPDIR/test_stdin.bed"
generate_bed 1000000 "$BENCH_FILE"

BENCH_STDIN=""
bench_one "file" "$PIO" "$BENCH_FILE"
file_ms=$RESULT_MS; file_kb=$RESULT_KB

BENCH_STDIN="$BENCH_FILE"
bench_one "stdin" "$PIO" "-"
stdin_ms=$RESULT_MS; stdin_kb=$RESULT_KB
BENCH_STDIN=""

printf "  File:  %s, %s peak RSS\n" "$(fmt_time $file_ms)" "$(fmt_size $file_kb)"
printf "  Stdin: %s, %s peak RSS\n" "$(fmt_time $stdin_ms)" "$(fmt_size $stdin_kb)"

echo ""

# ---------------------------------------------------------------------------
# Convert results CSV to README plot format and regenerate PNGs via plot_readme.gp
# ---------------------------------------------------------------------------

# benchmark_results.csv: 15 cols, time-ms and memory-kb interleaved per tool.
# benchmark_readme.csv:  15 cols, reads + 7 ms cols + 7 mb cols. Used by plot_readme.gp.
README_CSV="$SCRIPT_DIR/benchmark_readme.csv"
awk -F',' 'BEGIN{OFS=","}
NR==1 {
    print "reads,pio1_ms,pio8_ms,piolm_ms,sort1_ms,sort8_ms,bt_ms,bo_ms,pio1_mb,pio8_mb,piolm_mb,sort1_mb,sort8_mb,bt_mb,bo_mb"
    next
}
{
    printf "%s,%d,%d,%d,%d,%d,%d,%d,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f\n",
        $1, $2,$4,$6,$8,$10,$12,$14,
        $3/1024,$5/1024,$7/1024,$9/1024,$11/1024,$13/1024,$15/1024
}' "$CSV_FILE" > "$README_CSV"
echo "Wrote $README_CSV"

if command -v gnuplot &>/dev/null; then
    (cd "$SCRIPT_DIR" && gnuplot plot_readme.gp)
    echo "Plots regenerated in $SCRIPT_DIR (4 PNGs: time/memory × log/linear)"
else
    echo "gnuplot not found — skipping plot generation"
fi

echo ""
echo "Done."
