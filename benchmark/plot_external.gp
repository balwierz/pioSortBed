# Plot wall time, peak RSS, and total disk writes vs. input size
# for the --external-merge / --multi-pass / bedops / GNU sort
# comparison. Reads benchmark/bench_external_<BUDGET>_t<THREADS>.csv.
#
# Usage:
#   gnuplot -e "BUDGET='4G';THREADS=8" plot_external.gp
#   gnuplot -e "BUDGET='16G';THREADS=8" plot_external.gp

if (!exists("BUDGET"))  BUDGET  = '4G'
if (!exists("THREADS")) THREADS = 8
csv = sprintf("benchmark/bench_external_%s_t%d.csv", BUDGET, THREADS)

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
# CSV columns since the threads-aware bench:
# 1=tool 2=threads 3=reads 4=input_bytes 5=wall_s 6=peak_rss_kb
# 7=fs_writes_blocks_512 8=fs_reads_blocks_512 9=fs_writes_bytes 10=fs_reads_bytes
system(sprintf("awk -F, 'NR>1 && $1==\"pio-extmerge-zstd\"' %s > benchmark/.bx_pioe.tmp", csv))
system(sprintf("awk -F, 'NR>1 && $1==\"pio-multipass\"'     %s > benchmark/.bx_piom.tmp", csv))
system(sprintf("awk -F, 'NR>1 && $1==\"bedops-sort-bed\"'   %s > benchmark/.bx_beops.tmp", csv))
system(sprintf("awk -F, 'NR>1 && $1==\"gnu-sort\"'          %s > benchmark/.bx_gnu.tmp", csv))

pio_label_e = sprintf('pio --external-merge -t %d (zstd)', THREADS)
pio_label_m = sprintf('pio --multi-pass -t %d', THREADS)
gnu_label   = sprintf('GNU sort -S %s --parallel=%d', BUDGET, THREADS)
bops_label  = 'bedops sort-bed (1 core)'

# --- Wall time ---
set output sprintf('benchmark/bench_external_%s_t%d_time.png', BUDGET, THREADS)
set title sprintf("Wall time @ %s cap, %d threads (log-x, linear-y)", BUDGET, THREADS) font ',22'
set ylabel 'Wall time (s)' font ',20'
set yrange [0:*]
set format y '%g'
plot \
    'benchmark/.bx_pioe.tmp' using 3:5 with linespoints ls 1 title pio_label_e, \
    'benchmark/.bx_piom.tmp' using 3:5 with linespoints ls 2 title pio_label_m, \
    'benchmark/.bx_beops.tmp' using 3:5 with linespoints ls 3 title bops_label, \
    'benchmark/.bx_gnu.tmp'   using 3:5 with linespoints ls 4 title gnu_label

# --- Peak RSS ---
set output sprintf('benchmark/bench_external_%s_t%d_rss.png', BUDGET, THREADS)
set title sprintf("Peak resident set size @ %s cap, %d threads", BUDGET, THREADS) font ',22'
set ylabel 'Peak RSS (GiB)' font ',20'
plot \
    'benchmark/.bx_pioe.tmp' using 3:($6/1048576) with linespoints ls 1 title pio_label_e, \
    'benchmark/.bx_piom.tmp' using 3:($6/1048576) with linespoints ls 2 title pio_label_m, \
    'benchmark/.bx_beops.tmp' using 3:($6/1048576) with linespoints ls 3 title bops_label, \
    'benchmark/.bx_gnu.tmp'   using 3:($6/1048576) with linespoints ls 4 title gnu_label

# --- Total disk writes ---
set output sprintf('benchmark/bench_external_%s_t%d_writes.png', BUDGET, THREADS)
set title sprintf("Total bytes written to disk (output + temps) @ %s cap, %d threads", BUDGET, THREADS) font ',22'
set ylabel 'Disk writes (GiB)' font ',20'
plot \
    'benchmark/.bx_pioe.tmp' using 3:($9/1073741824) with linespoints ls 1 title pio_label_e, \
    'benchmark/.bx_piom.tmp' using 3:($9/1073741824) with linespoints ls 2 title pio_label_m, \
    'benchmark/.bx_beops.tmp' using 3:($9/1073741824) with linespoints ls 3 title bops_label, \
    'benchmark/.bx_gnu.tmp'   using 3:($9/1073741824) with linespoints ls 4 title gnu_label

# Cleanup temp files.
system("rm -f benchmark/.bx_*.tmp")
