#! /bin/bash

gnuplot <<EOF
#set size square
set grid
# set xrange[1:512]
# set yrange[-0.05:0.8]
#set xtics 0.1
set xlabel "time"
set ylabel "flow rate"
#set key left
plot "stats" u 1:3 w lines
set term pdf
set out 'tmp.pdf'
plot "stats" u 1:3 w lines
exit
EOF
