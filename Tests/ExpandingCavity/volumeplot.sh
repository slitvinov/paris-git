#!/bin/bash


gnuplot <<EOF
set term gif
set output "Volume growth.gif"
set xlabel " time "
set ylabel " cavity volume "
plot "stats" u 1:11 notitle
EOF