set datafile separator ','
# Treat 0 values as "missing" so configurations skipped at huge sizes (e.g.
# pioSortBed -t 8 and bedtools at 100M+ where they'd OOM on this fixture)
# don't drop to log(0) = -inf.
set datafile missing '0'
set logscale x 10
set format x '%.0s%c'
set xtics nomirror
set ytics nomirror
set key top left font ',13' spacing 1.3
set grid xtics ytics lt 0 lw 0.5 lc rgb '#dddddd'
set xrange [1e4:]

set style line 1 lc rgb '#e41a1c' lw 2.2 pt 7  ps 1.0 # pioSortBed 1t (red, circle)
set style line 2 lc rgb '#377eb8' lw 2.2 pt 5  ps 1.0 # pioSortBed 8t (blue, square)
set style line 3 lc rgb '#f781bf' lw 2.2 pt 13 ps 1.0 # pioSortBed low-mem (pink, diamond)
set style line 4 lc rgb '#4daf4a' lw 2.2 pt 9  ps 1.0 # GNU sort 1t (green, triangle up)
set style line 5 lc rgb '#ff7f00' lw 2.2 pt 11 ps 1.0 # GNU sort 8t (orange, triangle down)
set style line 6 lc rgb '#984ea3' lw 2.2 pt 3  ps 1.0 # bedtools (purple, plus)
set style line 7 lc rgb '#a65628' lw 2.2 pt 8  ps 1.0 # bedops (brown, star)

# --- Wall time plot (log scale) ---
set terminal pngcairo size 900,600 enhanced font 'Arial,11'
set output 'benchmark_time.png'
set title "Wall time" font ',14'
set xlabel 'Number of reads'
set ylabel 'Wall time (ms)'
set logscale y 10
set format y '%.2g'

plot 'benchmark_readme.csv' \
    using 1:2 with linespoints ls 1 title 'pioSortBed 1t', \
    '' using 1:3 with linespoints ls 2 title 'pioSortBed 8t', \
    '' using 1:4 with linespoints ls 3 title 'pioSortBed low-mem SSD', \
    '' using 1:5 with linespoints ls 4 title 'GNU sort 1t', \
    '' using 1:6 with linespoints ls 5 title 'GNU sort 8t', \
    '' using 1:7 with linespoints ls 6 title 'bedtools sort', \
    '' using 1:8 with linespoints ls 7 title 'bedops sort-bed'

# --- Peak memory plot (log scale) ---
set output 'benchmark_memory.png'
set title "Peak memory (RSS)" font ',14'
set ylabel 'Peak RSS (MB)'
set logscale y 10
set format y '%.2g'

plot 'benchmark_readme.csv' \
    using 1:9  with linespoints ls 1 title 'pioSortBed 1t', \
    '' using 1:10 with linespoints ls 2 title 'pioSortBed 8t', \
    '' using 1:11 with linespoints ls 3 title 'pioSortBed low-mem SSD', \
    '' using 1:12 with linespoints ls 4 title 'GNU sort 1t', \
    '' using 1:13 with linespoints ls 5 title 'GNU sort 8t', \
    '' using 1:14 with linespoints ls 6 title 'bedtools sort', \
    '' using 1:15 with linespoints ls 7 title 'bedops sort-bed'

# --- Wall time plot (linear scale) ---
set output 'benchmark_time_linear.png'
set title "Wall time (linear scale)" font ',14'
set ylabel 'Wall time (ms)'
unset logscale y
set format y '%.0f'
set yrange [0:]

plot 'benchmark_readme.csv' \
    using 1:2 with linespoints ls 1 title 'pioSortBed 1t', \
    '' using 1:3 with linespoints ls 2 title 'pioSortBed 8t', \
    '' using 1:4 with linespoints ls 3 title 'pioSortBed low-mem SSD', \
    '' using 1:5 with linespoints ls 4 title 'GNU sort 1t', \
    '' using 1:6 with linespoints ls 5 title 'GNU sort 8t', \
    '' using 1:7 with linespoints ls 6 title 'bedtools sort', \
    '' using 1:8 with linespoints ls 7 title 'bedops sort-bed'

# --- Peak memory plot (linear scale) ---
set output 'benchmark_memory_linear.png'
set title "Peak memory (RSS, linear scale)" font ',14'
set ylabel 'Peak RSS (MB)'
unset logscale y
set format y '%.0f'
set yrange [0:]

plot 'benchmark_readme.csv' \
    using 1:9  with linespoints ls 1 title 'pioSortBed 1t', \
    '' using 1:10 with linespoints ls 2 title 'pioSortBed 8t', \
    '' using 1:11 with linespoints ls 3 title 'pioSortBed low-mem SSD', \
    '' using 1:12 with linespoints ls 4 title 'GNU sort 1t', \
    '' using 1:13 with linespoints ls 5 title 'GNU sort 8t', \
    '' using 1:14 with linespoints ls 6 title 'bedtools sort', \
    '' using 1:15 with linespoints ls 7 title 'bedops sort-bed'
