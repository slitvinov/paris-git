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
set nolabel
set xrange[-0.1:1.1]
set yrange[-0.05:0.15]
plot "out/output_location000" w lp notitle, "out/output_location001" w lp notitle, 0.5*x*(1-x) w l notitle, 0 w l notitle
set term pdf
set out 'tmp.pdf'
plot "out/output_location000" w lp notitle, "out/output_location001" w lp notitle, 0.5*x*(1-x) w l notitle, 0 w l notitle
exit
EOF

# display tmp.pdf
