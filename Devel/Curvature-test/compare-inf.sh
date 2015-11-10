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

if [ $2 == 2 ]; then
gnuplot <<EOF
set log x
set log y
set xrange [2:110]
set xlabel "Grid points per Diameter"
set ylabel "Curvature error in $dim"
plot "gerris-$dim.txt" u (2*\$1):2 t "Gerris L2 one case" w lp,  "paris-$dim-$ndepth.tmp" u (2*\$1):2 t "ParisSim average L2", 4/(x*x) t "x^-2",  "paris-$dim-$ndepth.tmp" u (2*\$1):3 t "ParisSim max Linf" w lp
EOF
else
gnuplot <<EOF
set log x
set log y
set xrange [2:110]
set xlabel "Grid points per Diameter"
set ylabel "Curvature error in $dim"
plot "gerris-$dim.txt" u (2*\$1):2 t "Gerris L2 one case" w lp,  "paris-$dim-$ndepth.tmp" u (2*\$1):2 t "ParisSim average L2", 4/(x*x) t "x^-2",  "paris-$dim-$ndepth.tmp" u (2*\$1):3 t "ParisSim max Linf" w lp, 2/x
EOF
fi
