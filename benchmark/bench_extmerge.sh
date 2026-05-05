#!/bin/bash
# Compare --external-merge codecs on a large synthetic BED6 fixture.
# Reports wall time, peak RSS, and total temp-file bytes per codec.
#
# Usage:
#   bash benchmark/bench_extmerge.sh [num_reads] [budget] [tmpdir]
# Defaults: 100,000,000 reads, --max-mem 256M, /var/tmp
#
# Requires a WITH_BAM=1 build for rans0/rans1 to be available.

set -e
cd "$(dirname "$0")/.."

NUM_READS="${1:-100000000}"
BUDGET="${2:-256M}"
TMPDIR_BENCH="${3:-/var/tmp}"

PIO=./pioSortBed
BENCH_FILE="$TMPDIR_BENCH/extmerge_bench_${NUM_READS}.bed"
RUN_TMP="$TMPDIR_BENCH/extmerge_runs"
mkdir -p "$RUN_TMP"

if [[ ! -x "$PIO" ]]; then
    echo "Error: $PIO not built. Try 'make' or 'make WITH_BAM=1 HTSLIB=...'" >&2
    exit 1
fi

# Detect rANS support (only available in WITH_BAM builds).
HAS_RANS=0
if "$PIO" --merge-codec rans0 --external-merge /dev/null 2>&1 | grep -q "requires a WITH_BAM"; then
    HAS_RANS=0
else
    HAS_RANS=1
fi

CODECS=(raw lz4 zstd)
if (( HAS_RANS )); then CODECS+=(rans0 rans1); fi

# Generate fixture if missing.
if [[ ! -f "$BENCH_FILE" ]]; then
    echo "Generating $NUM_READS-record fixture at $BENCH_FILE ..."
    python3 -c "
import random, sys
random.seed(42)
n = int(sys.argv[1])
chroms = [f'chr{i}' for i in range(1, 23)]
out = sys.stdout
for i in range(n):
    c = random.choice(chroms)
    s = random.randint(1, 200_000_000)
    e = s + random.randint(50, 500)
    strand = '+' if random.random() < 0.5 else '-'
    out.write(f'{c}\t{s}\t{e}\tread{i}\t{random.randint(0,100)}\t{strand}\n')
" "$NUM_READS" > "$BENCH_FILE"
fi

INPUT_SIZE=$(stat -c '%s' "$BENCH_FILE")
echo "Input: $BENCH_FILE  ($(numfmt --to=iec --suffix=B "$INPUT_SIZE"), $NUM_READS records)"
echo "Budget: $BUDGET   tmpdir: $RUN_TMP"
echo "Hardware: $(grep 'model name' /proc/cpuinfo | head -1 | sed 's/.*: //')"
echo ""

CSV="benchmark/extmerge_${NUM_READS}_${BUDGET}.csv"
echo "codec,wall_s,peak_rss_kb,temp_bytes,num_runs,wall_per_byte_ns,ratio_vs_raw,ratio_vs_input" > "$CSV"

RAW_TEMP_BYTES=0  # set by raw run, used to compute ratio for others

printf "%-7s  %8s  %12s  %12s  %5s  %10s  %10s\n" \
    "codec" "wall_s" "peak_rss" "temp_bytes" "runs" "vs_raw" "vs_input"
printf "%s\n" "-------  --------  ------------  ------------  -----  ----------  ----------"

for codec in "${CODECS[@]}"; do
    # Clean stale runs.
    rm -f "$RUN_TMP"/piosort.run.*.tmp

    # Drop OS file-cache for this fixture (best-effort; needs root for full).
    # We use posix_fadvise via dd if=... iflag=nocache to evict fixture pages.
    # Skip if it fails — the comparison is still fair since all codecs run
    # against the same warmed cache state.
    sync
    dd if="$BENCH_FILE" of=/dev/null bs=1M iflag=nocache count=1 2>/dev/null || true

    LOG=$(mktemp)
    /usr/bin/time -f '%e %M' -o "$LOG.time" \
        "$PIO" --external-merge --merge-codec "$codec" \
              --max-mem "$BUDGET" --tmpdir "$RUN_TMP" -v \
              "$BENCH_FILE" > /dev/null 2> "$LOG"

    wall=$(awk '{print $1}' "$LOG.time")
    peak=$(awk '{print $2}' "$LOG.time")
    # Verbose log line looks like:
    # "Pass 1: 13 runs, temp_bytes=148472008 (...) codec=lz4, 0 s"
    runs=$(grep -oP 'Pass 1: \K\d+' "$LOG" || echo 0)
    temp=$(grep -oP 'temp_bytes=\K\d+' "$LOG" || echo 0)

    if [[ "$codec" == "raw" ]]; then RAW_TEMP_BYTES=$temp; fi

    if (( temp > 0 && RAW_TEMP_BYTES > 0 )); then
        vs_raw=$(awk -v t=$temp -v r=$RAW_TEMP_BYTES 'BEGIN{printf "%.3fx", t/r}')
    else
        vs_raw="-"
    fi
    vs_input=$(awk -v t=$temp -v i=$INPUT_SIZE 'BEGIN{printf "%.3fx", t/i}')
    wall_per_byte=$(awk -v w=$wall -v i=$INPUT_SIZE 'BEGIN{printf "%.2f", w*1e9/i}')

    printf "%-7s  %8s  %10s KB  %12s  %5s  %10s  %10s\n" \
        "$codec" "$wall" "$peak" "$temp" "$runs" "$vs_raw" "$vs_input"
    echo "$codec,$wall,$peak,$temp,$runs,$wall_per_byte,$vs_raw,$vs_input" >> "$CSV"

    rm -f "$LOG" "$LOG.time"
done

# Cleanup stale runs (the last codec's might still be there if we crashed).
rm -f "$RUN_TMP"/piosort.run.*.tmp
echo ""
echo "Wrote: $CSV"
