#!/bin/bash
#set -x

./clean.sh
mpirun -np 9 paris > tmpout  
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

awk '{print $1,$6}' ./out/element-00001.dat > results.dat

gnuplot <<EOF
set xlabel "time"
set ylabel "vel bubble"
set term pdf
set out "rising-bubble.pdf"
plot "results.dat" title "VOF2LPP" w l, "ref_LPP.dat" title "LPP" w l, "ref_DNS.dat" title "DNS" w l
set term png
set out "rising-bubble.png"
replot
EOF

precision=1e-3
compare results.dat ref_D2L.dat $precision 1 1

GREEN="\\033[1;32m"
NORMAL="\\033[0m"

