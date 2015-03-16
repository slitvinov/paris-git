#! /bin/bash

# if [ $# -lt 1 ] 
# then
#  echo "usage   : " $0 "time"
#  echo "example : " $0 "100"
#  exit
# fi

gnuplot <<EOF
set size square
set grid
set log xy
set xrange[1:512]
# set yrange[-0.05:0.15]
set xtics 4
set xlabel "number of cores/processes"
set ylabel "cell updates per second"
set key left
plot "speed.txt" using 1:2 w lp t col , "speed.txt" u 1:3 w lp t col , "speed.txt" u 1:4 w lp t col, "speed.txt" u 1:5 w lp t col,  "speed.txt" u 1:6 w lp t col, "speed.txt" u 1:7 w lp t col,  "speed.txt" u 1:8 w lp t col, "speed.txt" u 1:9 w lp t col,  "speed.txt" u 1:10 w lp t col,   x*1e6 t "linear" 
exit
EOF

# display tmp.pdf
