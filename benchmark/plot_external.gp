# Plot wall time, peak RSS, and total disk writes vs. input size
# for the --external-merge / --multi-pass / bedops / GNU sort
# comparison. Reads benchmark/bench_external.csv.
#
# Usage:
#   cd benchmark && gnuplot plot_external.gp

set datafile separator ','
set datafile missing 'NA'
set key top left font ',16' spacing 1.2
set xtics nomirror font ',16'
set ytics nomirror font ',16'
set grid xtics ytics lt 0 lw 0.5 lc rgb '#dddddd'

set logscale x 10
set format x '%.0s%c'

set style line 1 lc rgb '#e41a1c' lw 4 dt 1 pt 7  ps 2  # pio --external-merge
set style line 2 lc rgb '#c51b7d' lw 4 dt 1 pt 13 ps 2  # pio --multi-pass
set style line 3 lc rgb '#a65628' lw 4 dt 1 pt 8  ps 2  # bedops sort-bed
set style line 4 lc rgb '#4daf4a' lw 4 dt 1 pt 9  ps 2  # GNU sort

set terminal pngcairo size 1500,1000 enhanced font 'Arial,18'
set xlabel 'Number of reads' font ',20'
set xrange [1e7:1e9]

# Build per-tool data files via inline filter; gnuplot has no native
# group-by-column, so we emit four temp files first.
system("awk -F, 'NR>1 && $1==\"pio-extmerge-zstd\"' benchmark/bench_external.csv > benchmark/.bx_pioe.tmp")
system("awk -F, 'NR>1 && $1==\"pio-multipass\"'     benchmark/bench_external.csv > benchmark/.bx_piom.tmp")
system("awk -F, 'NR>1 && $1==\"bedops-sort-bed\"'   benchmark/bench_external.csv > benchmark/.bx_beops.tmp")
system("awk -F, 'NR>1 && $1==\"gnu-sort\"'          benchmark/bench_external.csv > benchmark/.bx_gnu.tmp")

# --- Wall time ---
set output 'benchmark/bench_external_time.png'
set title "Wall time @ 16 GiB cap (log-x, linear-y)" font ',22'
set ylabel 'Wall time (s)' font ',20'
set yrange [0:*]
set format y '%g'
plot \
    'benchmark/.bx_pioe.tmp' using 2:4 with linespoints ls 1 title 'pio --external-merge (zstd)', \
    'benchmark/.bx_piom.tmp' using 2:4 with linespoints ls 2 title 'pio --multi-pass', \
    'benchmark/.bx_beops.tmp' using 2:4 with linespoints ls 3 title 'bedops sort-bed', \
    'benchmark/.bx_gnu.tmp'   using 2:4 with linespoints ls 4 title 'GNU sort -S 16G'

# --- Peak RSS ---
set output 'benchmark/bench_external_rss.png'
set title "Peak resident set size @ 16 GiB cap" font ',22'
set ylabel 'Peak RSS (GiB)' font ',20'
plot \
    'benchmark/.bx_pioe.tmp' using 2:($5/1048576) with linespoints ls 1 title 'pio --external-merge (zstd)', \
    'benchmark/.bx_piom.tmp' using 2:($5/1048576) with linespoints ls 2 title 'pio --multi-pass', \
    'benchmark/.bx_beops.tmp' using 2:($5/1048576) with linespoints ls 3 title 'bedops sort-bed', \
    'benchmark/.bx_gnu.tmp'   using 2:($5/1048576) with linespoints ls 4 title 'GNU sort -S 16G'

# --- Total disk writes ---
set output 'benchmark/bench_external_writes.png'
set title "Total bytes written to disk (output + temps) @ 16 GiB cap" font ',22'
set ylabel 'Disk writes (GiB)' font ',20'
plot \
    'benchmark/.bx_pioe.tmp' using 2:($8/1073741824) with linespoints ls 1 title 'pio --external-merge (zstd)', \
    'benchmark/.bx_piom.tmp' using 2:($8/1073741824) with linespoints ls 2 title 'pio --multi-pass', \
    'benchmark/.bx_beops.tmp' using 2:($8/1073741824) with linespoints ls 3 title 'bedops sort-bed', \
    'benchmark/.bx_gnu.tmp'   using 2:($8/1073741824) with linespoints ls 4 title 'GNU sort -S 16G'

# Cleanup temp files.
system("rm -f benchmark/.bx_*.tmp")
