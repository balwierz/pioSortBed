# Plot region-query latency per (reader × region size), warm-cache, 100M-row data.
# Grouped bar chart: x-axis = region size, color = reader, y = median latency (log).

CSV = "bench_query.csv"
DAT = "/tmp/bench_query_plot.dat"

# Extract just the (reader, size, cache_state, median_us) tuples we want for the
# main plot. Filter to warm-cache, 100M-row data (regions.bed.gz / regions.lociss)
# — keep the inputs comparable. Emit one section per reader in plot order.
system "awk -F, 'NR>1 && $4==\"warm\" && ($2==\"regions.bed.gz\" || $2==\"regions.lociss\"){print $1,$3,$6}' " . CSV . " > " . DAT . ".raw"
system "sort -k1,1 -k2,2n " . DAT . ".raw > " . DAT

set datafile separator " "

# ----------------------------------------------------------------------
# Wall-time bar chart, log y. One bar per (reader, size). Each size
# group has 4 bars (tabix, polars, duckdb, loci).
# ----------------------------------------------------------------------
set terminal pngcairo size 1100,650 enhanced font 'Arial,11'
set output "bench_query_latency.png"

set title "Region-query latency: tabix vs polars / DuckDB / Loci on a 100M-record sorted BED\n(warm cache, median of 333 queries)" font ',12' enhanced
set xlabel "Region size"
set ylabel "Median latency (microseconds, log scale; lower is better)"
set logscale y 10
set yrange [50:300000]
set grid ytics lt 0 lw 0.5 lc rgb '#cccccc'
set style data histograms
set style histogram clustered gap 1.2
set style fill solid 0.9 border lc rgb '#444444'
set boxwidth 0.95
set key top left font ',10'
set xtics rotate by 0 nomirror
set ytics nomirror

# Three size buckets, four readers each. Map: each line in DAT is
# "<reader> <size> <median_us>"; readers are sorted alphabetically by awk
# pre-sort, so the order in DAT is duckdb / loci / polars / tabix repeating
# 3 times. We'll explicitly extract per-reader rows.
system "awk '$1==\"tabix\"  {print $2,$3}' " . DAT . " > /tmp/q_tabix.dat"
system "awk '$1==\"polars\" {print $2,$3}' " . DAT . " > /tmp/q_polars.dat"
system "awk '$1==\"duckdb\" {print $2,$3}' " . DAT . " > /tmp/q_duckdb.dat"
system "awk '$1==\"loci\"   {print $2,$3}' " . DAT . " > /tmp/q_loci.dat"

# Use a clustered histogram. The x-tic labels come from the tabix file.
plot newhistogram "" lt 1, \
     "/tmp/q_tabix.dat"  using 2:xtic(stringcolumn(1)=="1000"?"1 kbp":stringcolumn(1)=="100000"?"100 kbp":"10 Mbp") title "tabix"   lc rgb '#888888', \
     "/tmp/q_polars.dat" using 2 title "polars"  lc rgb '#1f77b4', \
     "/tmp/q_duckdb.dat" using 2 title "DuckDB"  lc rgb '#2ca02c', \
     "/tmp/q_loci.dat"   using 2 title "Loci"    lc rgb '#ff7f0e'

set output "bench_query_latency.svg"
set terminal svg size 1100,650 enhanced font 'Arial,11'
replot
