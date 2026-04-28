#!/usr/bin/env bash
# Test suite for pioSortBed.
# Compares output against reference implementations (LC_ALL=C sort) and
# checks specific behaviours (header passthrough, natural sort, gzip, etc.).
#
# Usage: bash test/test.sh [path/to/pioSortBed]
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIO="${1:-$REPO_DIR/pioSortBed}"
TMPDIR_BASE=$(mktemp -d /tmp/pioSortBed-test.XXXXXX)
trap 'rm -rf "$TMPDIR_BASE"' EXIT

PASS=0
FAIL=0

ok()  { echo "  PASS: $1"; PASS=$((PASS+1)); }
fail(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

check_eq() {
    local label="$1" got="$2" want="$3"
    if [ "$got" = "$want" ]; then ok "$label"; else
        fail "$label"
        echo "    expected: $(echo "$want" | awk 'NR<=5')"
        echo "    got:      $(echo "$got"  | awk 'NR<=5')"
    fi
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Generate a synthetic BED3 file: $1=nlines, $2=nchroms, $3=maxcoord
gen_bed3() {
    local n=$1 nchr=$2 max=$3
    awk -v n="$n" -v nc="$nchr" -v m="$max" 'BEGIN{
        srand(42)
        for(i=0;i<n;i++){
            c=int(rand()*nc)+1
            s=int(rand()*m)
            e=s+int(rand()*200)+1
            print "chr"c"\t"s"\t"e
        }
    }'
}

# Generate BED6 with strand
gen_bed6() {
    local n=$1 nchr=$2 max=$3
    awk -v n="$n" -v nc="$nchr" -v m="$max" 'BEGIN{
        srand(42)
        for(i=0;i<n;i++){
            c=int(rand()*nc)+1
            s=int(rand()*m)
            e=s+int(rand()*200)+1
            strand=(rand()>0.5)?"+":"-"
            print "chr"c"\t"s"\t"e"\tread"i"\t1\t"strand
        }
    }'
}

# Reference sort: LC_ALL=C sort -k1,1 -k2,2n
ref_sort()   { LC_ALL=C sort -k1,1 -k2,2n; }
ref_sort_b() { LC_ALL=C sort -k1,1 -k2,2n -k3,3n; }

# Canonical form: sort all columns so tie-breaking order doesn't matter.
# Both pioSortBed and GNU sort are non-stable for equal primary keys.
canon() { LC_ALL=C sort; }

# ---------------------------------------------------------------------------
echo "=== pioSortBed test suite ==="
echo "Binary: $PIO"
if [[ ! -x "$PIO" ]]; then
    echo "Error: binary not found or not executable: $PIO" >&2; exit 1
fi

# ---------------------------------------------------------------------------
echo ""
echo "--- Basic sort (file input, classic sort path) ---"

BED3="$TMPDIR_BASE/basic.bed"
gen_bed3 5000 5 1000000 > "$BED3"

want_s=$(ref_sort < "$BED3" | canon)
got=$("$PIO" "$BED3" 2>/dev/null | canon)
check_eq "sort-s file" "$got" "$want_s"

got=$("$PIO" --sort s "$BED3" 2>/dev/null | canon)
check_eq "sort-s explicit" "$got" "$want_s"

# ---------------------------------------------------------------------------
echo ""
echo "--- Sort by end (--sort b) ---"

want_b=$(ref_sort_b < "$BED3" | canon)
got=$("$PIO" --sort b "$BED3" 2>/dev/null | canon)
check_eq "sort-b file" "$got" "$want_b"

# ---------------------------------------------------------------------------
echo ""
echo "--- Stdin input ---"

got=$(cat "$BED3" | "$PIO" - 2>/dev/null | canon)
check_eq "stdin basic" "$got" "$want_s"

# ---------------------------------------------------------------------------
echo ""
echo "--- Gzip input ---"

GZ="$TMPDIR_BASE/basic.bed.gz"
gzip -c "$BED3" > "$GZ"
got=$("$PIO" "$GZ" 2>/dev/null | canon)
check_eq "gzip input" "$got" "$want_s"

# ---------------------------------------------------------------------------
echo ""
echo "--- Header line passthrough ---"

HDR="$TMPDIR_BASE/header.bed"
{
    echo "# comment line"
    echo "track name=test description=\"test track\""
    echo "browser position chr1:1-1000"
    cat "$BED3"
} > "$HDR"

got=$("$PIO" "$HDR" 2>/dev/null)
first=$(echo "$got" | awk 'NR<=3')
want_hdr=$(printf '# comment line\ntrack name=test description="test track"\nbrowser position chr1:1-1000')
check_eq "header passthrough lines" "$first" "$want_hdr"

# Body should be sorted
got_body=$(echo "$got" | tail -n +4 | canon)
want_body=$(ref_sort < "$BED3" | canon)
check_eq "header passthrough body sorted" "$got_body" "$want_body"

# ---------------------------------------------------------------------------
echo ""
echo "--- Natural sort (--natural-sort) ---"

NAT="$TMPDIR_BASE/natural.bed"
{
    echo "chr10	100	200"
    echo "chr2	100	200"
    echo "chr1	100	200"
    echo "chr20	100	200"
    echo "chrX	100	200"
    echo "chrY	100	200"
    echo "chr9	100	200"
} > "$NAT"

got=$("$PIO" --natural-sort "$NAT" 2>/dev/null | awk '{print $1}' | tr '\n' ' ' | xargs)
want="chr1 chr2 chr9 chr10 chr20 chrX chrY"
check_eq "natural-sort order" "$got" "$want"

# Default (lexicographic) sort should differ for chr2 vs chr10
got_lex=$("$PIO" "$NAT" 2>/dev/null | awk '{print $1}' | tr '\n' ' ' | xargs)
want_lex="chr1 chr10 chr2 chr20 chr9 chrX chrY"
check_eq "lexicographic sort order" "$got_lex" "$want_lex"

# ---------------------------------------------------------------------------
echo ""
echo "--- Sort by 5-prime end (--sort 5) ---"

BED6="$TMPDIR_BASE/stranded.bed"
gen_bed6 3000 3 500000 > "$BED6"

# For + strand: sort by start; for - strand: sort by end.
# We verify output is stable within each chromosome by checking 5' coords are non-decreasing.
got=$("$PIO" --sort 5 "$BED6" 2>/dev/null)
# For each chromosome, check 5' positions (beg for +, end for -) are non-decreasing
monotonic=$(echo "$got" | awk 'BEGIN{ok=1}{pos=($6=="-")?$3:$2; if($1==last&&pos<lpos){ok=0;exit}; last=$1;lpos=pos}END{print ok}')
[[ "$monotonic" == "1" ]] && ok "sort-5 non-decreasing within chr" || fail "sort-5 non-decreasing within chr"

# ---------------------------------------------------------------------------
echo ""
echo "--- Bucket sort path (--bucket-cutoff 0) ---"

got_bucket=$("$PIO" --bucket-cutoff 0 "$BED3" 2>/dev/null | canon)
check_eq "bucket sort path" "$got_bucket" "$want_s"

# ---------------------------------------------------------------------------
echo ""
echo "--- Low-mem-ssd mode ---"

got_lm=$("$PIO" --low-mem-ssd "$BED3" 2>/dev/null | canon)
check_eq "low-mem-ssd sort" "$got_lm" "$want_s"

# ---------------------------------------------------------------------------
echo ""
echo "--- Empty file ---"

EMPTY="$TMPDIR_BASE/empty.bed"
touch "$EMPTY"
got=$("$PIO" "$EMPTY" 2>/dev/null)
check_eq "empty file" "$got" ""

# ---------------------------------------------------------------------------
echo ""
echo "--- Single line ---"

SINGLE="$TMPDIR_BASE/single.bed"
echo -e "chr1\t100\t200" > "$SINGLE"
got=$("$PIO" "$SINGLE" 2>/dev/null)
check_eq "single line" "$got" "$(cat "$SINGLE")"

# ---------------------------------------------------------------------------
echo ""
echo "--- Multi-chromosome stress (many chromosomes, check ordering) ---"

MULTI="$TMPDIR_BASE/multi.bed"
# chromosomes not in sorted order to stress the sort
for i in 5 3 1 4 2; do
    for j in 300 100 200; do
        echo -e "chr${i}\t${j}\t$((j+50))"
    done
done > "$MULTI"

got=$("$PIO" "$MULTI" 2>/dev/null | canon)
want=$(ref_sort < "$MULTI" | canon)
check_eq "multi-chromosome ordering" "$got" "$want"

# ---------------------------------------------------------------------------
echo ""
echo "--- Collapse mode ---"

COLLAPSE="$TMPDIR_BASE/collapse.bed"
{
    echo "chr1	100	101	.	1.5	+"
    echo "chr1	100	101	.	2.5	+"
    echo "chr1	200	201	.	3.0	+"
    echo "chr2	50	51	.	4.0	+"
} > "$COLLAPSE"

got=$("$PIO" --collapse "$COLLAPSE" 2>/dev/null)
# position 100 should have weight 4.0 (1.5+2.5), position 200 → 3.0
line1=$(echo "$got" | head -1)
line2=$(echo "$got" | sed -n '2p')
line3=$(echo "$got" | sed -n '3p')
check_eq "collapse chr1:100 weight" "$(echo "$line1" | awk '{print $5}')" "4"
check_eq "collapse chr1:200 weight" "$(echo "$line2" | awk '{print $5}')" "3"
check_eq "collapse chr2:50 weight"  "$(echo "$line3" | awk '{print $5}')" "4"

# ---------------------------------------------------------------------------
echo ""
echo "--- Parallel paths (-t 8) ---"

# Larger fixture so the parallel parser actually splits; default-sort + sort-b + sort-5.
PAR="$TMPDIR_BASE/parallel.bed"
gen_bed6 200000 5 1000000 > "$PAR"

want_par=$(ref_sort   < "$PAR" | canon)
want_par_b=$(ref_sort_b < "$PAR" | canon)

got=$("$PIO" -t 8 "$PAR" 2>/dev/null | canon)
check_eq "-t 8 default" "$got" "$want_par"

got=$("$PIO" -t 8 --sort b "$PAR" 2>/dev/null | canon)
check_eq "-t 8 --sort b" "$got" "$want_par_b"

# --sort 5 with -t 8: monotonicity (multiset already covered by classic).
got_par5=$("$PIO" -t 8 --sort 5 "$PAR" 2>/dev/null)
mono_par5=$(echo "$got_par5" | awk 'BEGIN{ok=1}{pos=($6=="-")?$3:$2; if($1==last&&pos<lpos){ok=0;exit}; last=$1;lpos=pos}END{print ok}')
[[ "$mono_par5" == "1" ]] && ok "-t 8 --sort 5 monotonic within chr" || fail "-t 8 --sort 5 monotonic within chr"

# ---------------------------------------------------------------------------
echo ""
echo "--- Forced bucket-sort path with all sort modes ---"

got=$("$PIO" --bucket-cutoff 0 "$PAR" 2>/dev/null | canon)
check_eq "bucket --sort s" "$got" "$want_par"

got=$("$PIO" --bucket-cutoff 0 --sort b "$PAR" 2>/dev/null | canon)
check_eq "bucket --sort b" "$got" "$want_par_b"

got_buck5=$("$PIO" --bucket-cutoff 0 --sort 5 "$PAR" 2>/dev/null)
mono_buck5=$(echo "$got_buck5" | awk 'BEGIN{ok=1}{pos=($6=="-")?$3:$2; if($1==last&&pos<lpos){ok=0;exit}; last=$1;lpos=pos}END{print ok}')
[[ "$mono_buck5" == "1" ]] && ok "bucket --sort 5 monotonic within chr" || fail "bucket --sort 5 monotonic within chr"

# Bucket vs classic must agree on multiset (with --sort b they must be deterministic-equivalent).
got_b_classic=$("$PIO" --sort b "$PAR" 2>/dev/null | canon)
got_b_bucket=$("$PIO" --bucket-cutoff 0 --sort b "$PAR" 2>/dev/null | canon)
check_eq "bucket vs classic (--sort b)" "$got_b_classic" "$got_b_bucket"

# Parallel bucket vs serial bucket must agree on multiset.
got_pb=$("$PIO" -t 8 --bucket-cutoff 0 "$PAR" 2>/dev/null | canon)
got_sb=$("$PIO" -t 1 --bucket-cutoff 0 "$PAR" 2>/dev/null | canon)
check_eq "parallel bucket vs serial bucket" "$got_pb" "$got_sb"

# ---------------------------------------------------------------------------
echo ""
echo "--- RAL input format ---"

# RAL fields: id chr strand+context beg end weight ...
RAL="$TMPDIR_BASE/data.ral"
awk 'BEGIN{srand(42); for(i=0;i<2000;i++){
    c="chr"((i%5)+1); strand=(rand()>0.5)?"+":"-"
    beg=int(rand()*1000000); endp=beg+int(rand()*200)+1
    printf "id%d\t%s\t%s\t%d\t%d\t%g\n", i, c, strand, beg, endp, int(rand()*100)/10
}}' > "$RAL"

# RAL output preserves the input format. Reference: sort by chr (col 2), beg (col 4).
want_ral=$(LC_ALL=C sort -k2,2 -k4,4n "$RAL" | canon)
got=$("$PIO" --ral "$RAL" 2>/dev/null | canon)
check_eq "ral default" "$got" "$want_ral"

got=$("$PIO" --ral -t 8 "$RAL" 2>/dev/null | canon)
check_eq "ral -t 8" "$got" "$want_ral"

# ---------------------------------------------------------------------------
echo ""
echo "--- Version flag ---"

ver=$("$PIO" --version 2>&1)
[[ "$ver" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] && ok "version flag format" || fail "version flag format (got: $ver)"

# ---------------------------------------------------------------------------
echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
