set datafile separator ','
# Configurations skipped at huge sizes (e.g. pioSortBed -t 4 / -t 8 and
# bedtools at 100M+, where they'd OOM on this fixture) are recorded as 'NA'
# in the CSV. Tell gnuplot so those points don't render at all.
set datafile missing 'NA'
set xtics nomirror font ',18'
set ytics nomirror font ',18'
set key top left font ',18' spacing 1.2
set grid xtics ytics lt 0 lw 0.5 lc rgb '#dddddd'

# Same colour per tool family; thread count by line style. -t 1 solid, -t 4
# dashed (dt 2), -t 8 dotted (dt 3). pioSortBed classic and pioSortBed
# low-mem are different algorithms, different colours. bedtools / bedops
# are single-threaded by design.
set style line 1  lc rgb '#e41a1c' lw 4 dt 1 pt 7  ps 2 # pio 1t (red, circle, solid)
set style line 2  lc rgb '#e41a1c' lw 4 dt 2 pt 7  ps 2 # pio 4t (red, circle, dashed)
set style line 3  lc rgb '#e41a1c' lw 4 dt 3 pt 7  ps 2 # pio 8t (red, circle, dotted)
set style line 4  lc rgb '#c51b7d' lw 4 dt 1 pt 13 ps 2 # pio low-mem 1t (magenta, diamond, solid)
set style line 5  lc rgb '#c51b7d' lw 4 dt 2 pt 13 ps 2 # pio low-mem 4t (magenta, diamond, dashed)
set style line 6  lc rgb '#c51b7d' lw 4 dt 3 pt 13 ps 2 # pio low-mem 8t (magenta, diamond, dotted)
set style line 7  lc rgb '#4daf4a' lw 4 dt 1 pt 9  ps 2 # GNU sort 1t (green, triangle up, solid)
set style line 8  lc rgb '#4daf4a' lw 4 dt 2 pt 9  ps 2 # GNU sort 4t (green, triangle up, dashed)
set style line 9  lc rgb '#4daf4a' lw 4 dt 3 pt 9  ps 2 # GNU sort 8t (green, triangle up, dotted)
set style line 10 lc rgb '#984ea3' lw 4 dt 1 pt 3  ps 2 # bedtools sort (purple, plus)
set style line 11 lc rgb '#a65628' lw 4 dt 1 pt 8  ps 2 # bedops sort-bed (brown, star)

set terminal pngcairo size 1800,1200 enhanced font 'Arial,22'

# ============================================================================
# Log-log plots
# ============================================================================
set logscale x 10
set logscale y 10
set format x '%.0s%c'
set format y '%.2g'
set xrange [1e4:]
set yrange [*:*]
set xlabel 'Number of reads' font ',22'

# --- Wall time (log-log) ---
set output 'benchmark_time.png'
set title "Wall time (log-log)" font ',24'
set ylabel 'Wall time (ms)' font ',22'

plot 'benchmark_readme.csv' \
    using 1:2  with linespoints ls 1  title 'pioSortBed 1t', \
    '' using 1:3  with linespoints ls 2  title 'pioSortBed 4t', \
    '' using 1:4  with linespoints ls 3  title 'pioSortBed 8t', \
    '' using 1:5  with linespoints ls 4  title 'pioSortBed low-mem 1t', \
    '' using 1:6  with linespoints ls 5  title 'pioSortBed low-mem 4t', \
    '' using 1:7  with linespoints ls 6  title 'pioSortBed low-mem 8t', \
    '' using 1:8  with linespoints ls 7  title 'GNU sort 1t', \
    '' using 1:9  with linespoints ls 8  title 'GNU sort 4t', \
    '' using 1:10 with linespoints ls 9  title 'GNU sort 8t', \
    '' using 1:11 with linespoints ls 10 title 'bedtools sort', \
    '' using 1:12 with linespoints ls 11 title 'bedops sort-bed'

# --- Peak memory (log-log) ---
set output 'benchmark_memory.png'
set title "Peak memory (RSS, log-log)" font ',24'
set ylabel 'Peak RSS (MB)' font ',22'

plot 'benchmark_readme.csv' \
    using 1:13 with linespoints ls 1  title 'pioSortBed 1t', \
    '' using 1:14 with linespoints ls 2  title 'pioSortBed 4t', \
    '' using 1:15 with linespoints ls 3  title 'pioSortBed 8t', \
    '' using 1:16 with linespoints ls 4  title 'pioSortBed low-mem 1t', \
    '' using 1:17 with linespoints ls 5  title 'pioSortBed low-mem 4t', \
    '' using 1:18 with linespoints ls 6  title 'pioSortBed low-mem 8t', \
    '' using 1:19 with linespoints ls 7  title 'GNU sort 1t', \
    '' using 1:20 with linespoints ls 8  title 'GNU sort 4t', \
    '' using 1:21 with linespoints ls 9  title 'GNU sort 8t', \
    '' using 1:22 with linespoints ls 10 title 'bedtools sort', \
    '' using 1:23 with linespoints ls 11 title 'bedops sort-bed'

# ============================================================================
# Linear-linear plots (both axes linear)
# ============================================================================
unset logscale x
unset logscale y
set format x '%.0s%c'
set format y '%.0s%c'
set xrange [0:*]
set yrange [0:*]
set key top left

# --- Wall time (linear) ---
set output 'benchmark_time_linear.png'
set title "Wall time (linear)" font ',24'
set ylabel 'Wall time (ms)' font ',22'

plot 'benchmark_readme.csv' \
    using 1:2  with linespoints ls 1  title 'pioSortBed 1t', \
    '' using 1:3  with linespoints ls 2  title 'pioSortBed 4t', \
    '' using 1:4  with linespoints ls 3  title 'pioSortBed 8t', \
    '' using 1:5  with linespoints ls 4  title 'pioSortBed low-mem 1t', \
    '' using 1:6  with linespoints ls 5  title 'pioSortBed low-mem 4t', \
    '' using 1:7  with linespoints ls 6  title 'pioSortBed low-mem 8t', \
    '' using 1:8  with linespoints ls 7  title 'GNU sort 1t', \
    '' using 1:9  with linespoints ls 8  title 'GNU sort 4t', \
    '' using 1:10 with linespoints ls 9  title 'GNU sort 8t', \
    '' using 1:11 with linespoints ls 10 title 'bedtools sort', \
    '' using 1:12 with linespoints ls 11 title 'bedops sort-bed'

# --- Peak memory (linear) ---
set output 'benchmark_memory_linear.png'
set title "Peak memory (RSS, linear)" font ',24'
set ylabel 'Peak RSS (MB)' font ',22'

plot 'benchmark_readme.csv' \
    using 1:13 with linespoints ls 1  title 'pioSortBed 1t', \
    '' using 1:14 with linespoints ls 2  title 'pioSortBed 4t', \
    '' using 1:15 with linespoints ls 3  title 'pioSortBed 8t', \
    '' using 1:16 with linespoints ls 4  title 'pioSortBed low-mem 1t', \
    '' using 1:17 with linespoints ls 5  title 'pioSortBed low-mem 4t', \
    '' using 1:18 with linespoints ls 6  title 'pioSortBed low-mem 8t', \
    '' using 1:19 with linespoints ls 7  title 'GNU sort 1t', \
    '' using 1:20 with linespoints ls 8  title 'GNU sort 4t', \
    '' using 1:21 with linespoints ls 9  title 'GNU sort 8t', \
    '' using 1:22 with linespoints ls 10 title 'bedtools sort', \
    '' using 1:23 with linespoints ls 11 title 'bedops sort-bed'
