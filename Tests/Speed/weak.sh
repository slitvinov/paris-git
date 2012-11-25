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
#set yrange[50000:5e5]
set xtics 4
set key left
set xlabel "number of cores"
set ylabel "cell updates per second and per core (Z/np)"
set key left
plot "weak.txt" using (\$1):(\$2/\$1) with lp  t col, "weak.txt" u (\$1):(\$3/\$1) with lp  t col, 1.5e5*x**(-0.1) t "x**(-0.1)"
set term pdf
set out 'weak-scaling.pdf'
plot "weak.txt" using (\$1):(\$2/\$1) with lp  t col, "weak.txt" u (\$1):(\$3/\$1) with lp  t col, 1.5e5*x**(-0.1) t "x**(-0.1)"
exit
EOF

# display tmp.pdf
