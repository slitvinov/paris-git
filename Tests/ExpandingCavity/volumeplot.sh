#!/bin/bash


gnuplot <<EOF
set xlabel " time "
set ylabel " cavity volume "
plot "stats" u 1:11 notitle
EOF