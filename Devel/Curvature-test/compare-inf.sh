#! /bin/bash
#set -x

cp gerris-3D.txt gerris-3D.tmp
cp ../../Tests/Curvature/paris-$1.tmp . 
gnuplot <<EOF
call "plotinf.gp"
EOF