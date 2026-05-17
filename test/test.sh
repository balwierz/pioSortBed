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
echo "--- Low-mem-ssd mode ---"

got_lm=$("$PIO" --low-mem-ssd "$BED3" 2>/dev/null | canon)
check_eq "low-mem-ssd sort" "$got_lm" "$want_s"

# ---------------------------------------------------------------------------
echo ""
echo "--- External merge sort (--external-merge) ---"

# Generous budget — single run.
got_em=$("$PIO" --external-merge "$BED3" 2>/dev/null | canon)
check_eq "external-merge single run" "$got_em" "$want_s"

# Tight budget forces multiple runs + k-way merge.
got_em_multi=$("$PIO" --external-merge --max-mem 4K "$BED3" 2>/dev/null | canon)
check_eq "external-merge multi run" "$got_em_multi" "$want_s"

# ---------------------------------------------------------------------------
echo ""
echo "--- Multi-pass no-writes sort (--multi-pass) ---"

# Generous budget — single group.
got_mp=$("$PIO" --multi-pass "$BED3" 2>/dev/null | canon)
check_eq "multi-pass single group" "$got_mp" "$want_s"

# Tight budget forces multiple K-passes (one bucket = ~24 KB; 64 KB
# budget forces ~one chrom per pass on the BED3 fixture).
got_mp_multi=$("$PIO" --multi-pass --max-mem 64K "$BED3" 2>/dev/null | canon)
check_eq "multi-pass K-pass" "$got_mp_multi" "$want_s"

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
echo "--- --bgzip / --tabix (WITH_HTSLIB build only) ---"

# Skip these three tests on a non-htslib build. We detect by asking
# pioSortBed itself: it prints the --bgzip option in its help only on
# WITH_HTSLIB builds. If absent, skip cleanly.
if "$PIO" --help 2>&1 | grep -q -- '--bgzip'; then
    if ! command -v tabix &>/dev/null || ! command -v bgzip &>/dev/null; then
        echo "  SKIP: tabix or bgzip not installed (htslib CLI tools needed for tests)"
    else
        TMPBGZ=$(mktemp /tmp/pio_bgzip_XXXX.bed.gz)
        rm -f "$TMPBGZ" "$TMPBGZ.tbi"

        # 1. --bgzip round-trip: bgzipped output decompresses to the
        #    canonical sorted reference.
        wantbgz=$(ref_sort < "$PAR" | canon)
        "$PIO" -t 4 --bgzip -o "$TMPBGZ" "$PAR" 2>/dev/null
        gotbgz=$(bgzip -d -c "$TMPBGZ" 2>/dev/null | canon)
        check_eq "--bgzip round-trip vs plain sort" "$gotbgz" "$wantbgz"
        rm -f "$TMPBGZ"

        # 2. --bgzip --tabix index correctness: a region query returns
        #    the same rows as the awk-extracted reference region.
        "$PIO" -t 4 --bgzip --tabix -o "$TMPBGZ" "$PAR" 2>/dev/null
        if [[ -f "$TMPBGZ.tbi" ]]; then
            # Pick the first chromosome that exists in the test fixture
            # and a wide region that should overlap at least one record.
            firstchr=$(ref_sort < "$PAR" | head -1 | awk '{print $1}')
            wanttbx=$(ref_sort < "$PAR" | awk -F'\t' -v c="$firstchr" \
                '$1==c && $2<2000000000 && $3>0 {print $1"\t"$2"\t"$3}' | canon)
            gottbx=$(tabix "$TMPBGZ" "$firstchr:0-2000000000" 2>/dev/null \
                | awk -F'\t' '{print $1"\t"$2"\t"$3}' | canon)
            check_eq "--bgzip --tabix region query matches awk reference" \
                     "$gottbx" "$wanttbx"
        else
            fail "--bgzip --tabix did not produce a .tbi sidecar"
        fi
        rm -f "$TMPBGZ" "$TMPBGZ.tbi"

        # 3. --bgzip without -o must fail with exit code 1.
        if "$PIO" --bgzip "$PAR" > /dev/null 2>&1; then
            fail "--bgzip without -o should fail but succeeded"
        else
            ok "--bgzip without -o errors out (as documented)"
        fi
    fi
else
    echo "  SKIP: pioSortBed built without WITH_HTSLIB; --bgzip not available"
fi

# ---------------------------------------------------------------------------
echo ""
echo "--- Version flag ---"

ver=$("$PIO" --version 2>&1)
[[ "$ver" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] && ok "version flag format" || fail "version flag format (got: $ver)"

# ---------------------------------------------------------------------------
echo ""
echo "=== Results: $PASS passed, $FAIL failed ==="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
