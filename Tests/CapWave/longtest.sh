#!/bin/bash
#set -x

tmp=`mktemp -d 2>/dev/null || mktemp -d -t 'tmp'`

rm -fR out
cp inputlong input
mpirun -np 3 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

GREEN="\\033[1;32m"
NORMAL="\\033[0m"

ampini=$(head -n 1 interface.dat | awk '{printf("%10.9f",$2 - 1.50)}')

awk '{print $1, ($2-1.50)/'$ampini'*0.01}' interface.dat > $tmp/sim
awk '{print $2}' prosperetti > $tmp/theory
paste $tmp/sim $tmp/theory > comparison.dat

err=$(awk 'BEGIN{sum = 0}{ sum=sum+(($3-$2)*($3-$2))}END{ print sqrt(sum/NR)}' comparison.dat)

awk '{if ('$err' < 0.001) {print "\033[32;1m PASS\033[0m L2 relative error norm:" '$err'} else {print "\033[31;1m FAIL\033[0m L2 relative error norm:" '$err'}}' comparison.dat | tail -1

rm -rf $tmp

gnuplot <<EOF
set xlabel "time"
set ylabel "normalized amplitude a(t)/lambda"
p "./interface.dat" u (\$1*11.1366559937):((\$2-1.50)/$ampini*0.01) every 10 w p t "Simulation", "./prosperetti" u 1:2 w l t "Reference"
set term pngcairo
set out "capwave.png"
replot
EOF

if [ -d ../Testreport ] ; then
    mv capwave.png ../Testreport
fi

