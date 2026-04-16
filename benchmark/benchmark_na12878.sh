#!/usr/bin/env bash
set -euo pipefail

# Benchmark pioSortBed on real NA12878 (HG001) WGS reads aligned to GRCh38.
# Source: GIAB/NHGRI 300x Illumina HiSeq BAM, streamed via samtools.
# Default region: chr20 (~128M reads, ~5 GB BED, ~12 GB streamed from NCBI FTP).
# Usage: bash benchmark_na12878.sh [REGION]
#   REGION: samtools region string, default "chr20"
#   Examples: "chr1"  "chr20"  "chr1 chr2 chr3"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIO="$REPO_DIR/pioSortBed"
TIME_BIN=/usr/bin/time

BAM_URL="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam"

# Region(s) to extract — passed as separate samtools arguments
if [[ $# -ge 1 ]]; then
    REGIONS=("$@")
else
    REGIONS=("chr20")
fi

REGION_TAG="${REGIONS[*]// /_}"
BED_FILE="$SCRIPT_DIR/NA12878_GRCh38_${REGION_TAG}.bed"

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------

if [[ ! -x "$PIO" ]]; then
    echo "Error: pioSortBed binary not found at $PIO (run 'make' first)" >&2; exit 1
fi
for tool in samtools bedtools sort sort-bed; do
    command -v "$tool" &>/dev/null || echo "Warning: $tool not found, will be skipped"
done

HAS_SORT=0; HAS_BEDOPS=0
command -v sort     &>/dev/null && HAS_SORT=1
command -v sort-bed &>/dev/null && HAS_BEDOPS=1

# ---------------------------------------------------------------------------
# Download / prepare BED file
# ---------------------------------------------------------------------------

if [[ ! -f "$BED_FILE" ]]; then
    # Rough size estimates for chr20 at 300x (scale accordingly for other regions)
    echo "NA12878 WGS BED file not found: $(basename "$BED_FILE")"
    echo ""
    echo "Dataset : HG001 (NA12878) 300x Illumina WGS, GRCh38 — GIAB/NHGRI"
    echo "Source  : NCBI FTP (samtools HTTP streaming)"
    echo "Region  : ${REGIONS[*]}"
    echo ""
    echo "Estimated transfer : ~12 GB (chr20; proportionally more for larger regions)"
    echo "Estimated BED size : ~5 GB (chr20)"
    echo "Requirements       : samtools, bedtools, internet access"
    echo ""
    read -r -p "Proceed with streaming and conversion? [y/N] " confirm
    [[ "$confirm" =~ ^[Yy]$ ]] || { echo "Aborted."; exit 0; }

    echo ""
    echo "Streaming BAM → BED (samtools | bedtools bamtobed)..."
    echo "This may take 10–30 minutes depending on network speed."
    SECONDS=0
    samtools view -b -@ 4 "$BAM_URL" "${REGIONS[@]}" \
        | bedtools bamtobed -i stdin \
        > "$BED_FILE"
    echo "Done in ${SECONDS}s. $(wc -l < "$BED_FILE") reads → $BED_FILE"
fi

READ_COUNT=$(wc -l < "$BED_FILE")
FILE_SIZE=$(du -h "$BED_FILE" | cut -f1)
echo ""
echo "BED file : $BED_FILE"
echo "Reads    : $READ_COUNT"
echo "Size     : $FILE_SIZE"
echo ""

TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/pioSortBed-na12878.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

SORT_BUF="80%"

# ---------------------------------------------------------------------------
# Helpers (same pattern as benchmark.sh)
# ---------------------------------------------------------------------------

run_and_measure() {
    local time_out="$TMPDIR/time.tmp"
    "$TIME_BIN" -f '%e %M' -o "$time_out" "$@" \
        > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp"
    read -r secs kb < "$time_out"
    RESULT_MS=$(awk "BEGIN { printf \"%d\", $secs * 1000 }")
    RESULT_KB=$kb
}

run_shell_and_measure() {
    local time_out="$TMPDIR/time.tmp"
    "$TIME_BIN" -f '%e %M' -o "$time_out" bash -c "$1" \
        > "$TMPDIR/output.tmp" 2>"$TMPDIR/stderr.tmp"
    read -r secs kb < "$time_out"
    RESULT_MS=$(awk "BEGIN { printf \"%d\", $secs * 1000 }")
    RESULT_KB=$kb
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

# Warm page cache
cat "$BED_FILE" > /dev/null

# ---------------------------------------------------------------------------
# Run tools (all read from the file directly; cache warmed above)
# ---------------------------------------------------------------------------

echo "Running benchmarks..."
echo ""

printf "%-22s  %14s  %10s\n" "Tool" "Wall time" "Peak RSS"
printf "%-22s  %14s  %10s\n" "----" "---------" "--------"

# pioSortBed 1-thread
run_and_measure "$PIO" -t 1 "$BED_FILE"
pio1_ms=$RESULT_MS; pio1_kb=$RESULT_KB
printf "%-22s  %14s  %10s\n" "pioSortBed 1t" "$(fmt_time $pio1_ms)" "$(fmt_size $pio1_kb)"

# pioSortBed 8-thread (≥50M reads uses bucket sort, threads have no effect)
run_and_measure "$PIO" -t 8 "$BED_FILE"
pio8_ms=$RESULT_MS; pio8_kb=$RESULT_KB
printf "%-22s  %14s  %10s\n" "pioSortBed 8t" "$(fmt_time $pio8_ms)" "$(fmt_size $pio8_kb)"

# pioSortBed low-mem-ssd (two-pass, requires seekable file)
run_and_measure "$PIO" --low-mem-ssd "$BED_FILE"
piolm_ms=$RESULT_MS; piolm_kb=$RESULT_KB
printf "%-22s  %14s  %10s\n" "pioSortBed low-mem" "$(fmt_time $piolm_ms)" "$(fmt_size $piolm_kb)"

# GNU sort 1-thread
if (( HAS_SORT )); then
    run_shell_and_measure "LC_ALL=C sort -k1,1 -k2,2n --parallel=1 --buffer-size=$SORT_BUF '$BED_FILE'"
    sort1_ms=$RESULT_MS; sort1_kb=$RESULT_KB
    printf "%-22s  %14s  %10s\n" "GNU sort 1t" "$(fmt_time $sort1_ms)" "$(fmt_size $sort1_kb)"
fi

# GNU sort 8-thread
if (( HAS_SORT )); then
    run_shell_and_measure "LC_ALL=C sort -k1,1 -k2,2n --parallel=8 --buffer-size=$SORT_BUF '$BED_FILE'"
    sort8_ms=$RESULT_MS; sort8_kb=$RESULT_KB
    printf "%-22s  %14s  %10s\n" "GNU sort 8t" "$(fmt_time $sort8_ms)" "$(fmt_size $sort8_kb)"
fi

# bedops sort-bed
if (( HAS_BEDOPS )); then
    run_shell_and_measure "sort-bed '$BED_FILE'"
    bo_ms=$RESULT_MS; bo_kb=$RESULT_KB
    printf "%-22s  %14s  %10s\n" "bedops sort-bed" "$(fmt_time $bo_ms)" "$(fmt_size $bo_kb)"
fi

# bedtools: skip — memory usage at 100M+ reads exceeds practical limits
echo ""
echo "Note: bedtools sort skipped (peak RSS exceeds ~40 GB at this read count)"

echo ""
echo "Dataset: NA12878 HG001 GRCh38 300x WGS, region: ${REGIONS[*]}, reads: $READ_COUNT"
echo "Done."
