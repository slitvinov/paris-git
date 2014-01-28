#!/bin/bash
#set -x

gnuplot <<EOF &
set size square
set grid
set nolabel
set title "Droplet velocity in 3d lid-driven cavity" 
set xlabel "t (s)"
set ylabel "v (m/s)"
set xrange [0:0.004]
set yrange [-0.1:0.15]
plot './out_DNS/element-00001.dat' using 1:5 title "u,DNS",'./out_LPP/element-00001.dat' using 1:5 title "u,LPP", './out_DNS/element-00001.dat' using 1:6 title "v,DNS",'./out_LPP/element-00001.dat' using 1:6 title "v,LPP", './out_DNS/element-00001.dat' using 1:7 title "w,DNS",'./out_LPP/element-00001.dat' using 1:7 title "w,LPP"
pause 2
set term pdf
set out 'drop-vel.pdf'
plot './out_DNS/element-00001.dat' using 1:5 title "u,DNS",'./out_LPP/element-00001.dat' using 1:5 title "u,LPP", './out_DNS/element-00001.dat' using 1:6 title "v,DNS",'./out_LPP/element-00001.dat' using 1:6 title "v,LPP", './out_DNS/element-00001.dat' using 1:7 title "w,DNS",'./out_LPP/element-00001.dat' using 1:7 title "w,LPP"
exit
EOF

