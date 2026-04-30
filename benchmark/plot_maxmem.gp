set datafile separator ','
set terminal pngcairo size 1500,900 enhanced font 'Arial,20'
set output 'maxmem_sweep.png'

set title "pioSortBed --low-mem-ssd -t 4 @ 200M reads: --max-mem budget sweep" font ',22'
set xlabel "--max-mem budget" font ',20'
set ylabel "Wall time (s)" font ',20' textcolor rgb '#e41a1c'
set y2label "Peak RSS (GB)" font ',20' textcolor rgb '#377eb8'

set xtics nomirror font ',16'
set ytics nomirror font ',16' textcolor rgb '#e41a1c'
set y2tics font ',16' textcolor rgb '#377eb8'

set yrange [0:45]
set y2range [10:18]
set logscale x 2
set xrange [200e6:24e9]
set xtics ('256M' 268435456, '1G' 1073741824, '4G' 4294967296, '16G' 17179869184)
set grid xtics ytics lt 0 lw 0.5 lc rgb '#dddddd'
set key bottom right font ',16' spacing 1.3

# Uncapped reference values (no --max-mem flag): 20.16 s / 16084.4 MB
WALL_UNCAPPED_S = 20.16
RSS_UNCAPPED_GB = 16084.4 / 1024.0

# Skip header line; plot wall vs budget (left axis) + RSS vs budget (right axis)
# plus dashed reference lines for the uncapped baseline.
set datafile commentschars "#"
plot 'maxmem_sweep.csv' every ::1 using 2:($3/1000) with linespoints \
        lc rgb '#e41a1c' lw 4 pt 7 ps 2 axes x1y1 \
        title 'Wall time', \
     'maxmem_sweep.csv' every ::1 using 2:($5/1024) with linespoints \
        lc rgb '#377eb8' lw 4 pt 5 ps 2 axes x1y2 \
        title 'Peak RSS', \
     WALL_UNCAPPED_S with lines lc rgb '#e41a1c' lw 2 dt 2 axes x1y1 \
        title sprintf('Wall (uncapped, %.1f s)', WALL_UNCAPPED_S), \
     RSS_UNCAPPED_GB with lines lc rgb '#377eb8' lw 2 dt 2 axes x1y2 \
        title sprintf('RSS (uncapped, %.1f GB)', RSS_UNCAPPED_GB)
