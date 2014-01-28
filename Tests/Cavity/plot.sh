#!/bin/bash
#set -x

gnuplot <<EOF &
set size square
set grid
set nolabel
set title "Trajectory of droplet in 3d lid-driven cavity" 
set xlabel "x (m)"
set ylabel "y (m)"
set xrange [0:0.00032]
set yrange [0:0.00032]
plot './out_DNS/element-00001.dat' using 2:3 title "DNS",'./out_LPP/element-00001.dat' using 2:3 title "LPP"
pause 4
set term pdf
set out 'drop-traj.pdf'
plot './out_DNS/element-00001.dat' using 2:3 title "DNS",'./out_LPP/element-00001.dat' using 2:3 title "LPP"
exit
EOF

