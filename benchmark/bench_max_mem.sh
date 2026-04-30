#!/usr/bin/env bash
set -euo pipefail

# pio-lm -t 4 at 200M reads, sweeping --max-mem from tight (256M) to
# uncapped to map out the time/memory trade-off curve.

PIO=/home/piotr/Sources/pioSortBed/pioSortBed
TMPDIR=/var/tmp/pioSortBed-maxmem.$$
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr10" "chr11" "chr20" "chrX" "chrY")
N=200000000
FIXTURE="$TMPDIR/test_200M.bed"
echo "Generating 200M-row fixture (~8.6 GB)..." >&2
awk -v n="$N" -v nchr=${#CHROMS[@]} '
BEGIN {
    split("'"$(IFS=,; echo "${CHROMS[*]}")"'", chrs, ",")
    srand(42)
    for (i = 0; i < n; i++) {
        c = chrs[int(rand() * nchr) + 1]
        beg = int(rand() * 249000000)
        end = beg + int(rand() * 1000) + 1
        printf "%s\t%d\t%d\tread%d\t%d\t%s\n", c, beg, end, i, int(rand()*1000), (rand()>0.5 ? "+" : "-")
    }
}' > "$FIXTURE"
ls -lh "$FIXTURE" >&2
echo "" >&2

# (--max-mem value, label)
declare -a budgets=(
    "256M"
    "512M"
    "1G"
    "2G"
    "4G"
    "6G"
    "8G"
    "12G"
    "16G"
    "0"      # uncapped (no --max-mem flag)
)

OUT=/tmp/maxmem_results.csv
echo "max_mem,wall_ms,rss_kb,rss_mb" > "$OUT"

# Warm the page cache once
cat "$FIXTURE" > /dev/null

printf "%-10s  %12s  %10s\n" "max_mem" "wall_ms" "rss"
printf "%-10s  %12s  %10s\n" "-------" "-------" "---"

for budget in "${budgets[@]}"; do
    cat "$FIXTURE" > /dev/null  # keep cache warm
    if [[ "$budget" = "0" ]]; then
        cmd_args=( --low-mem-ssd -t 4 "$FIXTURE" )
        label="uncapped"
    else
        cmd_args=( --low-mem-ssd -t 4 --max-mem="$budget" "$FIXTURE" )
        label="$budget"
    fi

    /usr/bin/time -f '%e %M' -o "$TMPDIR/time.tmp" "$PIO" "${cmd_args[@]}" > /dev/null 2>"$TMPDIR/stderr.tmp"
    read -r secs kb < "$TMPDIR/time.tmp"
    ms=$(awk "BEGIN { printf \"%d\", $secs * 1000 }")
    mb=$(awk "BEGIN { printf \"%.1f\", $kb / 1024 }")
    printf "%-10s  %10s ms  %7.1f MB\n" "$label" "$ms" "$mb"
    echo "$label,$ms,$kb,$mb" >> "$OUT"
done

echo "" >&2
echo "Wrote $OUT" >&2
