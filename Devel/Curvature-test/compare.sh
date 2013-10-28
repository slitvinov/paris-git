#! /bin/bash
#set -x

cp gerris-3D.txt gerris-3D.tmp
cp ../../Tests/Curvature/paris-3D-?.tmp . 
gnuplot <<EOF
call "plot.gp"
EOF