#!/bin/bash
#set -x

if [ $# -lt 1 ]; then
    echo "missing arguments"
    echo usage $0 nx
    exit
fi

./run-2Dtest.sh 4 1 2e-4 $1 F 5e-4
./run-2Dtest.sh 4 1 2e-4 $1 T 5e-4


gnuplot <<EOF
set xlabel "time step"
set ylabel "flow rate"
#set key left
set log x
plot "flowrates-IMP-F.txt" u 1:3 w lp, "flowrates-IMP-T.txt" u 1:3 w lp
set term pdf
set out 'graph-nx-$1.pdf'
plot "flowrates-IMP-F.txt" u 1:3 w lp, "flowrates-IMP-T.txt" u 1:3 w lp
exit
EOF
