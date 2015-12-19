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
plot  4/(x*x) t "4 x^-2", "WithPopinet/popinet.csv" u (2*\$1):2 t "Popinet, 2009 Fig. 5" w lp, "paris-$dim-$ndepth.tmp" u (2*\$1):3 t "ParisSimulator" w lp,  "paris-$dim-$ndepth-nocent.tmp" u (2*\$1):3 t "ParisSim max Linf no centroids" w lp
EOF
else
gnuplot <<EOF
set log x
set log y
set grid
set xrange [1:110]
set xlabel "Grid points per Diameter"
set ylabel "Curvature error in $dim"
plot 2/(x*x) t "2/x",  "paris-$dim-$ndepth-nomixed.tmp" u (2*\$1):3 t "ParisSim max HF + centroids" w lp, "paris-$dim-$ndepth-nocent.tmp" u (2*\$1):3 t "ParisSim max HF + HF fit" w lp
EOF
fi

