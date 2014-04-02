#!/bin/bash


gnuplot <<EOF
set term gif
set output "Volume growth.gif"
set xlabel " time "
set ylabel " cavity volume "
set xrange [0:1]
phi(t) = pi*(0.2*0.2) + 0.05*4*t
plot "stats" u 1:11, phi(x)
set output "vectors_end.gif"
set yrange [0:1]
set xlabel " x "
set ylabel " y "
unset key
plot './out/UV-00000-20000.txt' with vectors
EOF
