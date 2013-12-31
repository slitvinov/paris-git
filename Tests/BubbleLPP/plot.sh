#!/bin/bash
#set -x

gnuplot <<EOF &
set size square
set grid
set nolabel
set title "Velocity of a single bubble" 
set xlabel "time step (dt=5e-6s)"
set ylabel "v (m/s)"
set yrange [0:0.1]
plot 'out_DNS/element-00001.dat' using 1:6 title "DNS", 'out_LPP/element-00001.dat' using 1:6 title "LPP"
pause 4
set term pdf
set out 'bubble_vel.pdf'
plot 'out_DNS/element-00001.dat' using 1:6 title "DNS", 'out_LPP/element-00001.dat' using 1:6 title "LPP"
exit
EOF

echo -e "$GREEN" "PASS & Check results in bubble_vel.pdf."  "$NORMAL"


