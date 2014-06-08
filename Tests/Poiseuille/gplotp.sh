#! /bin/bash
#set -x

if [ $# -lt 1 ] 
then
tail=`ls out/UV-* | awk 'BEGIN {FS="-"} {print $3}' | tail -1`
else
tail=$1.txt
fi

if  [ $# -lt 2 ] 
then
  rm -f 2plot.txt
  cat  out/UV-0000*-$tail >> 2plot.txt   
else
    cp out/UV-0000$2-$tail 2plot.txt
fi

if  [ $# -lt 2 ] 
then
  rm -f presplot.txt
  cat  out/P-0000*-$tail >> presplot.txt   
else
    cp out/P-0000$2-$tail presplot.txt
fi
gnuplot <<EOF
load "gnuplot.gp"
EOF

cp -p presplot.txt presplot.tmp
rm -f 2plot.txt presplot.txt
