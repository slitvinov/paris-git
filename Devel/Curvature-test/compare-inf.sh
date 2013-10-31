#! /bin/bash
#set -x


if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0  ndepth dim
    exit
fi

ndepth=$1
dim=$2'D'

cp gerris-3D.txt gerris-3D.tmp
cp ../../Tests/Curvature/paris-$dim-$ndepth.tmp . 
gnuplot <<EOF
set log x
set log y
set xlabel "Grid points per Radius"
set ylabel "Curvature error in $dim"
plot "gerris-$dim.txt" u 1:2 t "Gerris L2 one case",  "paris-$dim-$ndepth.tmp" u 1:2 t "ParisSim average L2", 1/(x*x) t "x^-2",  "paris-$dim-$ndepth.tmp" u 1:3 t "ParisSim max Linf", 1/x t "x^-1"
set term pdf
set out "curvature-$dim-$ndepth.pdf"
plot "gerris-$dim.txt" u 1:2 t "Gerris L2 one case",  "paris-$dim-$ndepth.tmp" u 1:2 t "ParisSim average L2", 1/(x*x) t "x^-2",  "paris-$dim-$ndepth.tmp" u 1:3 t "ParisSim max Linf", 1/x t "x^-1"
EOF