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
if [ -f ../../Tests/Curvature/paris-$dim-$ndepth.tmp ] ; then
    cp ../../Tests/Curvature/paris-$dim-$ndepth.tmp .
else
    echo "Warning no fresh paris test found"
fi
cp paris-$dim-$ndepth.tmp 1.tmp
egrep -v nan 1.tmp > paris-$dim-$ndepth.tmp
rm -f 1.tmp



if [ $2 == 2 ]; then
gnuplot <<EOF
set log x
set log y
set grid
set xrange [1:100]
set yrange [1e-3:10]
set xlabel "Grid points per Diameter"
set ylabel "Curvature relative L_inf error in $dim"
plot  4/(x*x) t "4 x^-2", "WithPopinet/popinet.csv" u (2*\$1):2 t "Popinet, 2009 Fig. 5 max norm" w lp,  "paris-$dim-$ndepth.tmp" u (2*\$1):3 t "ParisSim max max" w lp
set term pngcairo size 600, 600
set out "curvature$dim.png"
set size square
replot
EOF
else
gnuplot <<EOF
set log x
set log y
set grid
set xrange [1:110]
set xlabel "Grid points per Diameter"
set ylabel "Curvature Linf error norm in $dim"
plot 2/(x*x) t "2/(x*x)", "bsk-$dim-log.txt" u 1:4 t "Basilisk-$dim" w lp, "paris-$dim-$ndepth.tmp" u (2*\$1):3 t "ParisSim current" w lp, "paris-$dim-$ndepth-old.txt" u (2*\$1):3 t "ParisSim old method" w lp
set term pngcairo size 600, 600
set out "curvature$dim.png"
set size square
replot
EOF
fi

