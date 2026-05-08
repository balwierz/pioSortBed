#!/usr/bin/env bash
set -euo pipefail

# Benchmark all pioSortBed sort modes on real NA12878 (HG001) WGS data,
# compare against the previous version (v3.2.0 in pioSortBed.bak/) and
# external tools (GNU sort, bedops sort-bed, bedtools sort).
#
# Default region: all100M (100M reads sampled across all standard chroms).
# Streaming modes (--external-merge, --multi-pass) only support the
# default coordinate sort, so --sort=b/5 are tested only on classic path.
#
# Usage: bash benchmark_na12878.sh [REGION]
#   REGION: matches NA12878_GRCh38_${REGION}.bed (default: "all100M")
# Override the v3.2.0 binary with PIO_OLD=/path/to/binary

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIO="$REPO_DIR/pioSortBed"
PIO_OLD="${PIO_OLD:-$REPO_DIR/../pioSortBed.bak/pioSortBed}"
TIME_BIN=/usr/bin/time

REGION="${1:-all100M}"
BED_FILE="$SCRIPT_DIR/NA12878_GRCh38_${REGION}.bed"
MAX_MEM="${MAX_MEM:-4G}"

[[ -x "$PIO" ]]     || { echo "Error: $PIO not found (run 'make pioSortBed')" >&2; exit 1; }
[[ -f "$BED_FILE" ]] || { echo "Error: $BED_FILE not found" >&2; exit 1; }

PIO_OLD_AVAIL=0
if [[ -x "$PIO_OLD" ]]; then
    PIO_OLD_VER=$("$PIO_OLD" --version 2>&1 | head -1)
    PIO_OLD_AVAIL=1
fi
PIO_NEW_VER=$("$PIO" --version 2>&1 | head -1)

HAS_SORT=0; HAS_BEDOPS=0; HAS_BEDTOOLS=0
command -v sort     &>/dev/null && HAS_SORT=1
command -v sort-bed &>/dev/null && HAS_BEDOPS=1
command -v bedtools &>/dev/null && HAS_BEDTOOLS=1

READ_COUNT=$(wc -l < "$BED_FILE")
# Follow symlinks (-L) so we get the real file size, not the symlink's 71 bytes.
FILE_SIZE_H=$(du -hL "$BED_FILE" | cut -f1)
FILE_SIZE_B=$(stat -Lc %s "$BED_FILE")

echo "==============================================================="
echo "BED file       : $BED_FILE"
echo "Reads          : $READ_COUNT"
echo "Size           : $FILE_SIZE_H ($FILE_SIZE_B bytes)"
echo "Current binary : $PIO ($PIO_NEW_VER)"
if (( PIO_OLD_AVAIL )); then
    echo "Previous bin   : $PIO_OLD ($PIO_OLD_VER)"
else
    echo "Previous bin   : NOT FOUND ($PIO_OLD) — version comparison skipped"
fi
echo "max-mem cap    : $MAX_MEM (for --low-mem-ssd / --external-merge / --multi-pass)"
echo "==============================================================="
echo ""

TMPDIR_BENCH=$(mktemp -d "${TMPDIR:-/tmp}/pioSortBed-na12878.XXXXXX")
trap 'rm -rf "$TMPDIR_BENCH"' EXIT

CSV="$SCRIPT_DIR/benchmark_na12878_${REGION}.csv"
echo "tool,version,wall_s,peak_rss_kb" > "$CSV"

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

run_one() {
    # run_one <label> <version> <cmd...>
    local label="$1" version="$2"; shift 2
    local time_out="$TMPDIR_BENCH/time.tmp"
    "$TIME_BIN" -f '%e %M' -o "$time_out" "$@" \
        > "$TMPDIR_BENCH/output.tmp" 2>"$TMPDIR_BENCH/stderr.tmp" || {
        echo "  $label: FAILED (exit $?)" >&2
        cat "$TMPDIR_BENCH/stderr.tmp" >&2 | tail -10
        echo "$label,$version,FAIL,FAIL" >> "$CSV"
        return
    }
    read -r secs kb < "$time_out"
    echo "$label,$version,$secs,$kb" >> "$CSV"
    printf "  %-32s  %10ss  %10s\n" "$label" "$secs" "$(fmt_size $kb)"
}

run_shell() {
    local label="$1" version="$2"; shift 2
    local time_out="$TMPDIR_BENCH/time.tmp"
    "$TIME_BIN" -f '%e %M' -o "$time_out" bash -c "$1" \
        > "$TMPDIR_BENCH/output.tmp" 2>"$TMPDIR_BENCH/stderr.tmp" || {
        echo "  $label: FAILED" >&2
        echo "$label,$version,FAIL,FAIL" >> "$CSV"
        return
    }
    read -r secs kb < "$time_out"
    echo "$label,$version,$secs,$kb" >> "$CSV"
    printf "  %-32s  %10ss  %10s\n" "$label" "$secs" "$(fmt_size $kb)"
}

fmt_size() {
    awk -v kb="$1" 'BEGIN {
        if (kb >= 1048576) printf "%.1f GB", kb/1048576
        else if (kb >= 1024) printf "%.1f MB", kb/1024
        else printf "%d kB", kb
    }'
}

# Warm page cache so cold-disk effects don't dominate the first run.
echo "Warming page cache..."
cat "$BED_FILE" > /dev/null
echo ""

# -----------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------

printf "%-34s  %11s  %10s\n" "Tool / mode" "Wall time" "Peak RSS"
printf "%-34s  %11s  %10s\n" "------------------------------" "----------" "--------"

echo "--- Current pioSortBed ($PIO_NEW_VER), --sort=s (default) ---"
run_one "pio-classic-1t"       "$PIO_NEW_VER" "$PIO" -t 1 "$BED_FILE"
run_one "pio-classic-8t"       "$PIO_NEW_VER" "$PIO" -t 8 "$BED_FILE"
run_one "pio-low-mem-ssd-1t"   "$PIO_NEW_VER" "$PIO" --low-mem-ssd -t 1 --max-mem "$MAX_MEM" "$BED_FILE"
run_one "pio-low-mem-ssd-8t"   "$PIO_NEW_VER" "$PIO" --low-mem-ssd -t 8 --max-mem "$MAX_MEM" "$BED_FILE"
run_one "pio-external-merge-8t" "$PIO_NEW_VER" "$PIO" --external-merge -t 8 --max-mem "$MAX_MEM" --tmpdir "$TMPDIR_BENCH" "$BED_FILE"
run_one "pio-multi-pass-8t"    "$PIO_NEW_VER" "$PIO" --multi-pass -t 8 --max-mem "$MAX_MEM" "$BED_FILE"
echo ""

echo "--- Current pioSortBed ($PIO_NEW_VER), other sort orders (8t classic) ---"
run_one "pio-sort=b-8t"        "$PIO_NEW_VER" "$PIO" --sort b -t 8 "$BED_FILE"
run_one "pio-sort=5-8t"        "$PIO_NEW_VER" "$PIO" --sort 5 -t 8 "$BED_FILE"
echo ""

if (( PIO_OLD_AVAIL )); then
    echo "--- Previous pioSortBed ($PIO_OLD_VER) — version comparison ---"
    run_one "pio-OLD-classic-1t"     "$PIO_OLD_VER" "$PIO_OLD" -t 1 "$BED_FILE"
    run_one "pio-OLD-classic-8t"     "$PIO_OLD_VER" "$PIO_OLD" -t 8 "$BED_FILE"
    run_one "pio-OLD-low-mem-ssd-1t" "$PIO_OLD_VER" "$PIO_OLD" --low-mem-ssd "$BED_FILE"
    echo ""
fi

echo "--- External tools ---"
if (( HAS_SORT )); then
    run_shell "gnu-sort-1t"  "$(sort --version | head -1)" "LC_ALL=C sort -k1,1 -k2,2n --parallel=1 --buffer-size=80% '$BED_FILE'"
    run_shell "gnu-sort-8t"  "$(sort --version | head -1)" "LC_ALL=C sort -k1,1 -k2,2n --parallel=8 --buffer-size=80% '$BED_FILE'"
fi
if (( HAS_BEDOPS )); then
    # sort-bed prints version under --help, not --version; use a non-fatal pipeline
    BO_VER=$(sort-bed --help 2>&1 | awk -F'version: *' '/version:/ {print $2; exit}' || true)
    run_one "bedops-sort-bed" "${BO_VER:-bedops}" sort-bed "$BED_FILE"
fi
if (( HAS_BEDTOOLS )); then
    BT_VER=$(bedtools --version 2>&1 | head -1 || true)
    run_one "bedtools-sort"   "$BT_VER" bedtools sort -i "$BED_FILE"
fi

echo ""
echo "==============================================================="
echo "Done. CSV: $CSV"
echo "==============================================================="
