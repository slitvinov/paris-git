#! /bin/bash

if [ $# -lt 1 ] 
then
 echo "usage   : " $0 "itimestep  [processor rank]"
 echo "example : " $0 "01 0"
 exit
fi

if  [ $# -lt 2 ] 
then
  rm -f 2plot.txt
  cat  out/UV-00000-000$1.txt >> 2plot.txt   
  cat  out/UV-00001-000$1.txt >> 2plot.txt   
  cat  out/UV-00002-000$1.txt >> 2plot.txt   
  cat  out/UV-00003-000$1.txt >> 2plot.txt   
else
    cp out/UV-0000$2-000$1.txt 2plot.txt
fi
gnuplot <<EOF
set size square
set nolabel
#set xrange[-0.5:0.5]
#set yrange[-0.5:0.5]
plot "2plot.txt" w vec notitle
EOF