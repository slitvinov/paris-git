#!/bin/bash
#set -x

if [ $# -lt 1 ]; then
    echo "missing arguments"
    echo usage $0 nx
    exit
fi

#./run-2dtest.sh 6  2e-4 $1 F 5e-4
#./run-2dtest.sh 6  2e-4 $1 T 5e-4

./run-2dtest.sh 4  4e-4 $1 F 5e-4
./run-2dtest.sh 4  4e-4 $1 T 5e-4


cat > gnuplot.tmp <<EOF
plot "flowrates-IMP-F.txt" u 1:3 w lp t "explicit", "flowrates-IMP-T.txt" u 1:3 w lp t "implicit"
EOF

gnuplot <<EOF
set xlabel "time step"
set ylabel "flow rate"
#set key left
set log x
load "gnuplot.tmp"
set term pdf
set out 'graph-nx-$1.pdf'
load "gnuplot.tmp"
exit
EOF
