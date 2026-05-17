#!/usr/bin/env bash
# Benchmark pioSortBed BAM input vs `samtools sort` on the same BAM.
#
# Reports wall time, peak RSS, and total bytes written (output BAM + temp
# files). The bytes-written column is the SSD-wear-relevant metric; it
# captures samtools sort's spill behaviour when input ≫ -m budget.
#
# Usage: bench_bam.sh INPUT.bam [THREADS]
#   THREADS: passed to both tools. Default 8.
#
# Output: bench_bam.csv (one row per (tool × config))
#
set -euo pipefail
cd "$(dirname "$0")/.."

INPUT="${1:?'usage: bench_bam.sh INPUT.bam [THREADS]'}"
THREADS="${2:-8}"
[[ -f "$INPUT" ]] || { echo "Missing input: $INPUT" >&2; exit 1; }

PIO="$(pwd)/pioSortBed"
[[ -x "$PIO" ]] || { echo "Build first with: make WITH_BAM=1" >&2; exit 1; }

WORK=$(mktemp -d /var/tmp/bambench.XXXX)
trap 'rm -rf "$WORK"' EXIT

CSV="benchmark/bench_bam.csv"
echo "tool,threads,reads,input_bytes,wall_s,peak_rss_kb,fs_writes_blocks_512,fs_writes_bytes" > "$CSV"

READS=$(samtools view -c "$INPUT" 2>/dev/null)
INPUT_BYTES=$(stat -c '%s' "$INPUT")
echo "Input: $INPUT  ($(numfmt --to=iec --suffix=B $INPUT_BYTES), $READS reads)"
echo "Threads: $THREADS"
echo ""

drop_fcache() { sync; dd if="$1" of=/dev/null bs=1M iflag=nocache count=1 2>/dev/null || true; }

run() {
    local label="$1" outpath="$2"
    shift 2
    rm -f "$outpath"
    drop_fcache "$INPUT"
    local terr=$(mktemp)
    /usr/bin/time -v -o "$terr" "$@" 2>/dev/null

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
    printf "  %-32s wall=%-7s s peak=%-9s KB fs_writes=%-8s MB (out: %s)\n" \
        "$label" "$wall" "$peak" "$wb_mb" "$(stat -c '%s' "$outpath" 2>/dev/null | numfmt --to=iec --suffix=B)"
    echo "$label,$THREADS,$READS,$INPUT_BYTES,$wall,$peak,$wb,$wb_bytes" >> "$CSV"
}

# 1. samtools sort, default settings (768 MiB per thread).
run "samtools-sort-default" "$WORK/sam_default.bam" \
    samtools sort -@ "$THREADS" -o "$WORK/sam_default.bam" "$INPUT"

# 2. samtools sort with -m 4G (so it does NOT spill on smaller inputs; matches
#    pio's in-RAM design more directly).
run "samtools-sort-m4G" "$WORK/sam_m4G.bam" \
    samtools sort -@ "$THREADS" -m 4G -o "$WORK/sam_m4G.bam" "$INPUT"

# 3. pioSortBed (in-RAM, no spilling).
run "pioSortBed-t${THREADS}" "$WORK/pio.bam" \
    "$PIO" -t "$THREADS" -o "$WORK/pio.bam" "$INPUT"

echo ""
echo "Wrote: $CSV"

# Correctness: samtools view text output of each sorted file should be
# identical (modulo @PG lines which differ by tool). Sort the records to
# canonicalise (same coord ties resolve differently between tools).
echo ""
echo "Correctness check (samtools view | LC_ALL=C sort | md5):"
for f in "$WORK/sam_default.bam" "$WORK/sam_m4G.bam" "$WORK/pio.bam"; do
    md=$(samtools view "$f" 2>/dev/null | LC_ALL=C sort -k3,3 -k4,4n -k1,1 | md5sum | cut -d' ' -f1)
    printf "  %-32s %s\n" "$(basename "$f")" "$md"
done
