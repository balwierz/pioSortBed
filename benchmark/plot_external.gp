# Plot wall time, peak RSS, and total disk writes vs. input size
# for the --external-merge / --multi-pass / bedops / GNU sort
# comparison. Reads benchmark/bench_external_<BUDGET>.csv.
#
# Usage:
#   cd benchmark && gnuplot -e "BUDGET='16G'" plot_external.gp
#   cd benchmark && gnuplot -e "BUDGET='4G'" plot_external.gp

if (!exists("BUDGET")) BUDGET = '16G'
csv = sprintf("benchmark/bench_external_%s.csv", BUDGET)

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
system(sprintf("awk -F, 'NR>1 && $1==\"pio-extmerge-zstd\"' %s > benchmark/.bx_pioe.tmp", csv))
system(sprintf("awk -F, 'NR>1 && $1==\"pio-multipass\"'     %s > benchmark/.bx_piom.tmp", csv))
system(sprintf("awk -F, 'NR>1 && $1==\"bedops-sort-bed\"'   %s > benchmark/.bx_beops.tmp", csv))
system(sprintf("awk -F, 'NR>1 && $1==\"gnu-sort\"'          %s > benchmark/.bx_gnu.tmp", csv))

gnu_label = sprintf('GNU sort -S %s', BUDGET)

# --- Wall time ---
set output sprintf('benchmark/bench_external_%s_time.png', BUDGET)
set title sprintf("Wall time @ %s cap (log-x, linear-y)", BUDGET) font ',22'
set ylabel 'Wall time (s)' font ',20'
set yrange [0:*]
set format y '%g'
plot \
    'benchmark/.bx_pioe.tmp' using 2:4 with linespoints ls 1 title 'pio --external-merge (zstd)', \
    'benchmark/.bx_piom.tmp' using 2:4 with linespoints ls 2 title 'pio --multi-pass', \
    'benchmark/.bx_beops.tmp' using 2:4 with linespoints ls 3 title 'bedops sort-bed', \
    'benchmark/.bx_gnu.tmp'   using 2:4 with linespoints ls 4 title gnu_label

# --- Peak RSS ---
set output sprintf('benchmark/bench_external_%s_rss.png', BUDGET)
set title sprintf("Peak resident set size @ %s cap", BUDGET) font ',22'
set ylabel 'Peak RSS (GiB)' font ',20'
plot \
    'benchmark/.bx_pioe.tmp' using 2:($5/1048576) with linespoints ls 1 title 'pio --external-merge (zstd)', \
    'benchmark/.bx_piom.tmp' using 2:($5/1048576) with linespoints ls 2 title 'pio --multi-pass', \
    'benchmark/.bx_beops.tmp' using 2:($5/1048576) with linespoints ls 3 title 'bedops sort-bed', \
    'benchmark/.bx_gnu.tmp'   using 2:($5/1048576) with linespoints ls 4 title gnu_label

# --- Total disk writes ---
set output sprintf('benchmark/bench_external_%s_writes.png', BUDGET)
set title sprintf("Total bytes written to disk (output + temps) @ %s cap", BUDGET) font ',22'
set ylabel 'Disk writes (GiB)' font ',20'
plot \
    'benchmark/.bx_pioe.tmp' using 2:($8/1073741824) with linespoints ls 1 title 'pio --external-merge (zstd)', \
    'benchmark/.bx_piom.tmp' using 2:($8/1073741824) with linespoints ls 2 title 'pio --multi-pass', \
    'benchmark/.bx_beops.tmp' using 2:($8/1073741824) with linespoints ls 3 title 'bedops sort-bed', \
    'benchmark/.bx_gnu.tmp'   using 2:($8/1073741824) with linespoints ls 4 title gnu_label

# Cleanup temp files.
system("rm -f benchmark/.bx_*.tmp")
