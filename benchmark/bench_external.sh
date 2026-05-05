#!/bin/bash
# Benchmark pioSortBed --external-merge and --multi-pass against bedops
# sort-bed and GNU sort, with a fixed memory cap (default 16 GiB), as a
# function of input size. Reports wall time, peak RSS, and total bytes
# written to disk (via /proc/PID/io write_bytes).
#
# Usage:
#   bash benchmark/bench_external.sh "<sizes>" <budget> <tmpdir>
# Defaults: sizes "10000000 50000000 100000000 200000000 500000000"
#           budget 16G
#           tmpdir /var/tmp

set -e
cd "$(dirname "$0")/.."

SIZES="${1:-10000000 50000000 100000000 200000000 500000000}"
BUDGET="${2:-16G}"
TMPDIR_BENCH="${3:-/var/tmp}"

PIO=./pioSortBed
RUN_TMP="$TMPDIR_BENCH/extbench_runs"
mkdir -p "$RUN_TMP"
export TMPDIR="$RUN_TMP"
export LC_ALL=C

if [[ ! -x "$PIO" ]]; then
    echo "Error: $PIO not built. make first." >&2
    exit 1
fi

CSV="benchmark/bench_external.csv"
echo "tool,reads,input_bytes,wall_s,peak_rss_kb,fs_writes_blocks_512,fs_reads_blocks_512,fs_writes_bytes,fs_reads_bytes" > "$CSV"

# Run a command, capture wall, peak RSS, and File system inputs/outputs
# via /usr/bin/time -v. ru_oublock / ru_inblock count block-level I/O
# (1 block = 512 bytes); this is what the process *caused* to be queued
# at the block layer, including writes that may still be dirty in the
# page cache when the process exits. Reliable for short jobs unlike
# /proc/PID/io's write_bytes which only counts pages actually flushed
# during the process's lifetime.
# Result: prints "wall peak fs_writes fs_reads" on stdout.
run_with_io() {
    local out_file="$1"
    shift
    local err_file=$(mktemp)

    /usr/bin/time -v -o "$err_file" \
        "$@" > "$out_file" 2>/dev/null

    local wall=$(awk -F': ' '/Elapsed \(wall clock\)/ {print $NF}' "$err_file")
    # Convert h:mm:ss / m:ss.ss to seconds.
    wall=$(echo "$wall" | awk -F: '{
        if (NF==3) print $1*3600 + $2*60 + $3;
        else if (NF==2) print $1*60 + $2;
        else print $1
    }')
    local peak=$(awk -F': ' '/Maximum resident set size/ {print $NF}' "$err_file")
    local wb=$(awk -F': ' '/File system outputs/ {print $NF}' "$err_file")
    local rb=$(awk -F': ' '/File system inputs/  {print $NF}' "$err_file")
    rm -f "$err_file"
    echo "$wall $peak ${wb:-0} ${rb:-0}"
}

# Drop file-cache for a path (best-effort, no root needed).
drop_fcache() {
    sync
    dd if="$1" of=/dev/null bs=1M iflag=nocache count=1 2>/dev/null || true
}

for n in $SIZES; do
    BENCH_FILE="$TMPDIR_BENCH/bench_${n}.bed"
    OUT_FILE="$TMPDIR_BENCH/bench_${n}.out"

    if [[ ! -f "$BENCH_FILE" ]]; then
        echo ">>> generating $n-record BED6 fixture at $BENCH_FILE ..."
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
" "$n" > "$BENCH_FILE"
    fi
    INPUT_BYTES=$(stat -c '%s' "$BENCH_FILE")
    echo ""
    echo "=== reads=$n input=$(numfmt --to=iec --suffix=B $INPUT_BYTES) budget=$BUDGET ==="
    printf "%-22s  %8s  %10s  %14s  %14s\n" \
        "tool" "wall_s" "peak_RSS" "fs_writes_MB" "fs_reads_MB"
    printf "%s\n" "----------------------  --------  ----------  --------------  --------------"

    run_one() {
        local label="$1" csv_label="$2"; shift 2
        drop_fcache "$BENCH_FILE"
        rm -f "$RUN_TMP"/piosort.run.*.tmp 2>/dev/null
        local res wall peak wb rb
        res=$(run_with_io "$OUT_FILE" "$@")
        read wall peak wb rb <<<"$res"
        local wb_bytes=$((wb * 512))
        local rb_bytes=$((rb * 512))
        local wb_mb=$(awk -v b=$wb_bytes 'BEGIN{printf "%.1f", b/1048576}')
        local rb_mb=$(awk -v b=$rb_bytes 'BEGIN{printf "%.1f", b/1048576}')
        printf "%-22s  %8.2f  %8s KB  %14s  %14s\n" "$label" "$wall" "$peak" "$wb_mb" "$rb_mb"
        echo "$csv_label,$n,$INPUT_BYTES,$wall,$peak,$wb,$rb,$wb_bytes,$rb_bytes" >> "$CSV"
    }

    run_one "pio --external-merge" "pio-extmerge-zstd" \
        "$PIO" --external-merge --merge-codec zstd \
               --max-mem "$BUDGET" --tmpdir "$RUN_TMP" "$BENCH_FILE"

    run_one "pio --multi-pass" "pio-multipass" \
        "$PIO" --multi-pass --max-mem "$BUDGET" "$BENCH_FILE"

    run_one "bedops sort-bed" "bedops-sort-bed" \
        sort-bed --max-mem "$BUDGET" --tmpdir "$RUN_TMP" "$BENCH_FILE"

    run_one "GNU sort -S $BUDGET" "gnu-sort" \
        sort -k1,1 -k2,2n -S "$BUDGET" -T "$RUN_TMP" "$BENCH_FILE"

    rm -f "$OUT_FILE" "$RUN_TMP"/piosort.run.*.tmp
done

echo ""
echo "Wrote: $CSV"
echo ""
echo "Note: fs_writes_MB is the volume of block-layer writes the process"
echo "      caused (getrusage ru_oublock × 512 B). For SSD wear that's"
echo "      the right metric; pages still dirty in cache at process exit"
echo "      *do* count here (unlike /proc/PID/io write_bytes)."
