#! /bin/bash

# if [ $# -lt 1 ] 
# then
#  echo "usage   : " $0 "time"
#  echo "example : " $0 "100"
#  exit
# fi

gnuplot <<EOF
set grid
# set xrange[-0.1:1.1]
# set yrange[-0.05:0.15]
set xtics 4
set xlabel "number of cores/processes"
set ylabel "cell updates per second"
plot "speed.txt" using 1:2 w lp title "AMD 6136 16-core 'celsius'" , "speed.txt" u 1:3 w lp title "AMD 8378 32-core 'cooley'"
set term pdf
set out 'tmp.pdf'
plot "speed.txt" using 1:2 w lp title "AMD 6136" , "speed.txt" u 1:3 w lp title "AMD 8378"
exit
EOF

# display tmp.pdf
