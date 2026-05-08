# Plot v3.7.0 vs v3.5.0 streaming-mode performance.
# Bench_external sweep, 16 GiB budget, 8 threads, all tools.

CSV_NEW = "bench_external_16G_t8.csv"
CSV_OLD = "bench_external_16G_t8.v3.5.0_baseline.csv"

# Pre-filter into per-tool files so gnuplot's linespoints can connect
system "awk -F, 'NR>1 && $1==\"pio-extmerge-zstd\"{print $3,$5}' " . CSV_NEW . " > /tmp/pio-em-new.dat"
system "awk -F, 'NR>1 && $1==\"pio-multipass\"    {print $3,$5}' " . CSV_NEW . " > /tmp/pio-mp-new.dat"
system "awk -F, 'NR>1 && $1==\"bedops-sort-bed\"  {print $3,$5}' " . CSV_NEW . " > /tmp/bedops-new.dat"
system "awk -F, 'NR>1 && $1==\"gnu-sort\"         {print $3,$5}' " . CSV_NEW . " > /tmp/gnu-new.dat"
system "awk -F, 'NR>1 && $1==\"pio-extmerge-zstd\"{print $3,$5}' " . CSV_OLD . " > /tmp/pio-em-old.dat"
system "awk -F, 'NR>1 && $1==\"pio-multipass\"    {print $3,$5}' " . CSV_OLD . " > /tmp/pio-mp-old.dat"

set datafile separator " "

set terminal pngcairo size 1100,650 enhanced font 'Arial,11'
set output "bench_external_v37_vs_v35_time.png"

set title "External-merge & multi-pass — wall time vs input size\nv3.7.0 (current, solid) vs v3.5.0 (baseline, dotted) — 16 GiB budget, 8 threads" font ',12' enhanced
set xlabel "Reads"
set ylabel "Wall time (s, log scale)"
set logscale x 10
set logscale y 10
set format x "%.0s%c"
set grid xtics ytics lt 0 lw 0.5 lc rgb '#cccccc'
set key top left font ',10'

# v3.7.0 lines (solid)
# v3.5.0 lines (dashed, lighter)

plot "/tmp/pio-em-new.dat" using 1:2 with linespoints lc rgb '#1f77b4' lw 2.2 pt 7  ps 1.0 title "extmerge v3.7.0", \
     "/tmp/pio-mp-new.dat" using 1:2 with linespoints lc rgb '#2ca02c' lw 2.2 pt 7  ps 1.0 title "multipass v3.7.0", \
     "/tmp/pio-em-old.dat" using 1:2 with linespoints lc rgb '#1f77b4' lw 1.4 dt 3 pt 6 ps 1.0 title "extmerge v3.5.0", \
     "/tmp/pio-mp-old.dat" using 1:2 with linespoints lc rgb '#2ca02c' lw 1.4 dt 3 pt 6 ps 1.0 title "multipass v3.5.0", \
     "/tmp/bedops-new.dat" using 1:2 with linespoints lc rgb '#888888' lw 1.4 pt 5 ps 0.8 title "bedops sort-bed", \
     "/tmp/gnu-new.dat"    using 1:2 with linespoints lc rgb '#ff7f0e' lw 1.4 pt 5 ps 0.8 title "GNU sort -8t"

set output "bench_external_v37_vs_v35_time.svg"
set terminal svg size 1100,650 enhanced font 'Arial,11'
replot
