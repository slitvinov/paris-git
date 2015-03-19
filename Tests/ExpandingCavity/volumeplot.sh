#!/bin/bash

gnuplot <<EOF
set term gif
set output "Volume growth.gif"
set xlabel " time "
set ylabel " cavity volume "
phi(t) = pi*(0.2*0.2) + 0.05*4*t
plot "stats" u 1:11, phi(x)
EOF
