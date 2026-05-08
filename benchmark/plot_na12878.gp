# Plot NA12878 100M-read benchmark: wall time + peak RSS, sorted by speed.
# Bars are color-coded by category (v3.7.0 / v3.2.0 / external) using one
# plot command per category — gives us a proper gnuplot key for free.

CSV = "benchmark_na12878_all100M.csv"
DATAFILE = "/tmp/na12878_sorted.dat"

# Tab-delimited so labels with spaces stay one field. Columns:
#   1: x position (1..N), 2: wall_s, 3: peak_GB, 4: label, 5: category
system "awk -F, -v OFS='\\t' 'NR>1{ \
    cat = ($1 ~ /^pio-OLD/) ? 2 : (($1 ~ /^pio-/) ? 1 : 3); \
    label = $1; \
    sub(/^pio-OLD-/, \"v3.2.0 \", label); \
    sub(/^pio-/, \"v3.7.0 \", label); \
    print $3, $4/1048576, label, cat \
}' " . CSV . " | sort -t'\t' -k1,1n | awk -F'\\t' -v OFS='\\t' '{print NR, $0}' > " . DATAFILE

# ----------------------------------------------------------------------
# Common style
# ----------------------------------------------------------------------
set datafile separator "\t"
set boxwidth 0.7
set style fill solid 0.9 border lc rgb '#444444'
set grid ytics lt 0 lw 0.5 lc rgb '#cccccc'
set xtics rotate by -45 font ',9'
set key top left font ',10' box lc rgb '#cccccc'

C1 = '#1f77b4'  # v3.7.0 (current)
C2 = '#ff7f0e'  # v3.2.0 (baseline)
C3 = '#888888'  # external tools

# ----------------------------------------------------------------------
# 1. Wall time
# ----------------------------------------------------------------------
set terminal pngcairo size 1100,650 enhanced font 'Arial,11'
set output "benchmark_na12878_wall.png"

set title "NA12878 GRCh38 WGS — 100 M reads, 6.6 GB BED — wall time" font ',13'
set ylabel "Wall time (s, lower is better)"
set yrange [0:200]
set ytics 25
set xrange [0:16]

plot DATAFILE using 1:(0):xtic(stringcolumn(4)) with boxes lc rgb '#ffffff' notitle, \
     DATAFILE using ($5==1 ? $1 : 1/0):2 with boxes lc rgb C1 title "v3.7.0 (current)", \
     DATAFILE using ($5==2 ? $1 : 1/0):2 with boxes lc rgb C2 title "v3.2.0 (baseline)", \
     DATAFILE using ($5==3 ? $1 : 1/0):2 with boxes lc rgb C3 title "external tools", \
     DATAFILE using 1:($2+3):(sprintf("%.1fs", $2)) with labels font ',9' textcolor rgb '#333333' notitle

set output "benchmark_na12878_wall.svg"
set terminal svg size 1100,650 enhanced font 'Arial,11'
replot

# ----------------------------------------------------------------------
# 2. Peak RSS (same x-order)
# ----------------------------------------------------------------------
set terminal pngcairo size 1100,650 enhanced font 'Arial,11'
set output "benchmark_na12878_rss.png"

set title "NA12878 GRCh38 WGS — 100 M reads, 6.6 GB BED — peak RSS" font ',13'
set ylabel "Peak RSS (GB, lower is better)"
set yrange [0:55]
set ytics 5
set key top left

plot DATAFILE using 1:(0):xtic(stringcolumn(4)) with boxes lc rgb '#ffffff' notitle, \
     DATAFILE using ($5==1 ? $1 : 1/0):3 with boxes lc rgb C1 title "v3.7.0 (current)", \
     DATAFILE using ($5==2 ? $1 : 1/0):3 with boxes lc rgb C2 title "v3.2.0 (baseline)", \
     DATAFILE using ($5==3 ? $1 : 1/0):3 with boxes lc rgb C3 title "external tools", \
     DATAFILE using 1:($3+1.2):(sprintf("%.1f", $3)) with labels font ',9' textcolor rgb '#333333' notitle

set output "benchmark_na12878_rss.svg"
set terminal svg size 1100,650 enhanced font 'Arial,11'
replot
