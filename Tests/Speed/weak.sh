#! /bin/bash

file=$1
file2=$2

data_type=raw
if [ $data_type == raw ]; then

 if [ $# -lt 1 ] 
 then
  echo "usage   : $0 file"
  exit
 fi

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
plot "$file" using (\$1):(\$2) with lp  t col, "$file" u (\$1):(\$3) with lp  t col,  "$file" u (\$1):(\$4) with lp  t col, 1.5e5*x**(-0.2) t "x**(-0.2)"
set term pdf
set out 'weak-scaling.pdf'
plot "$file" using (\$1):(\$2) with lp  t col, "$file" u (\$1):(\$3) with lp  t col,  "$file" u (\$1):(\$4) with lp  t col, 1.5e5*x**(-0.2) t "x**(-0.2)"
exit
EOF
else

 if [ $# -lt 2 ] 
 then
  echo "usage   : $0 file1 file2"
  exit
 fi

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
plot "$file" using (\$1):(\$2/\$1) with lp  t col, "$file" u (\$1):(\$3/\$1) with lp  t col,  "$file2" u (\$1):(\$2/\$1) with lp  t col, 1.5e5*x**(-0.2) t "x**(-0.2)"
set term pdf
set out 'weak-scaling.pdf'
plot "$file" using (\$1):(\$2/\$1) with lp  t col, "$file" u (\$1):(\$3/\$1) with lp  t col,  "$file2" u (\$1):(\$2/\$1) with lp  t col, 1.5e5*x**(-0.2) t "x**(-0.2)"
exit
EOF
fi
# display tmp.pdf
