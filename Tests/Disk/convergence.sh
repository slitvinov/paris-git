#!/bin/bash
#set -x


awk  'BEGIN { reference=0.11269487E-01 } { print $1 " " ($3-reference)}' < flowrates-IMP-T.txt > 2plot.txt

cat > 2plot.gp <<EOF
plot "2plot.txt" w lp t "implicit method error", x
EOF
gnuplot <<EOF
set xlabel "dt"
set ylabel "error"
set log x
set log y
load "2plot.gp"
set term pdf
set out 'convergence.pdf'
load "2plot.gp"
exit
EOF
