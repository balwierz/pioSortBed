#!/usr/bin/env bash
# Pipeline-time benchmark: total wall + total bytes written for three
# end-to-end recipes that all produce a "sorted, indexed, queryable" BED.
#
# Recipes:
#   1) Canonical 3-tool: LC_ALL=C sort | bgzip > x.bed.gz ; tabix x.bed.gz
#   2) pio integrated:   pio --low-mem-ssd -t N --bgzip --tabix -o x.bed.gz
#   3) pio LociSSD:      pio --low-mem-ssd -t N --lociss-output x.lociss
#
# Each is run on the NA12878 100M-read BED. Output: pipeline.csv with
# (recipe, wall_s, peak_rss_kb, fs_writes_bytes).
#
set -euo pipefail
cd "$(dirname "$0")/.."

INPUT="${1:-benchmark/NA12878_GRCh38_all100M.bed}"
[[ -f "$INPUT" ]] || { echo "Missing input: $INPUT" >&2; exit 1; }

WORK=$(mktemp -d /var/tmp/pipeline.XXXX)
trap 'rm -rf "$WORK"' EXIT
THREADS=${THREADS:-8}

CSV="benchmark/bench_pipeline.csv"
echo "recipe,wall_s,peak_rss_kb,fs_writes_blocks_512,fs_writes_bytes" > "$CSV"

# Drop file cache (best-effort) so the input read isn't free
drop_fcache() { sync; dd if="$1" of=/dev/null bs=1M iflag=nocache count=1 2>/dev/null || true; }

# Run an inner command under /usr/bin/time -v, parse out wall/RSS/writes
# (single combined value across the whole pipeline — bash -c is parented).
run_recipe() {
    local label="$1" cmd="$2"
    local terr=$(mktemp)
    drop_fcache "$INPUT"
    echo "=== $label ==="
    /usr/bin/time -v -o "$terr" bash -c "$cmd"
    local wall peak wb
    wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {print $NF}' "$terr")
    wall=$(echo "$wall" | awk -F: '{
        if (NF==3) print $1*3600 + $2*60 + $3;
        else if (NF==2) print $1*60 + $2;
        else print $1
    }')
    peak=$(awk -F': ' '/Maximum resident set size/ {print $NF}' "$terr")
    wb=$(awk -F': ' '/File system outputs/ {print $NF}' "$terr")
    rm -f "$terr"
    local wb_bytes=$((wb * 512))
    local wb_mb
    wb_mb=$(awk -v b="$wb_bytes" 'BEGIN{printf "%.1f", b/1048576}')
    printf "  wall=%s s, peak_rss=%s KB, fs_writes=%s MB (%s blocks)\n" \
        "$wall" "$peak" "$wb_mb" "$wb"
    echo "$label,$wall,$peak,$wb,$wb_bytes" >> "$CSV"
}

run_recipe "canonical-sort-bgzip-tabix" \
    "LC_ALL=C sort -k1,1 -k2,2n --parallel=$THREADS -S 16G \"$INPUT\" \
     | bgzip --threads=$THREADS > \"$WORK/out.bed.gz\" && \
     tabix \"$WORK/out.bed.gz\""

run_recipe "pio-bgzip-tabix" \
    "./pioSortBed --low-mem-ssd -t $THREADS --bgzip --tabix \
     -o \"$WORK/pio.bed.gz\" \"$INPUT\""

run_recipe "pio-lociss" \
    "./pioSortBed --low-mem-ssd -t $THREADS \
     --lociss-output \"$WORK/pio.lociss\" \"$INPUT\""

run_recipe "pio-lociss-with-index" \
    "./pioSortBed --low-mem-ssd -t $THREADS \
     --lociss-output \"$WORK/pio_idx.lociss\" --lociss-index \"$INPUT\""

echo ""
echo "Wrote: $CSV"
