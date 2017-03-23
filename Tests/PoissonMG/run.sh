#!/bin/bash
#set -x

rm -rf tmpout out convergenceMG_history.txt

cp inputtemplate input
mpirun -np 2 paris > tmpout
res=$(tail -1 convergenceMG_history.txt | awk '{print $4} ')
cpu=$(tail -1 convergenceMG_history.txt | awk '{printf "%3.2g" , $3} ')
err=$(echo `awk ' /error/ { res = $2 } END { print res } ' < tmpout`)

#echo " "
echo "cpu =" $cpu 

awk '{if ('$res' < 0.00000001 && '$err' < 0.1) {print "\033[32;1m PASS\033[0m Residual '$res' Error '$err' "} 
         else {print "\033[31;1m FAIL\033[0m Residual '$res' Error '$err' "}}' tmpout | tail -1

gnuplot <<EOF 
set size square
set grid
set nolabel
set log y
set xlabel 'iteration'
set ylabel 'Error'
set y2label 'CPU time'
set key outside below
set pointsize 0.5
set y2tics
set ytics nomirror
plot "convergenceMG_history.txt" u 2:4 t 'residual' w lp pt 7, "convergenceMG_history.txt" u 2:3 t 'CPU time' ax x1y2 w p pt 7
pause 3
set term pdf
set out 'convergence.pdf'
replot
exit
EOF

rm -rf input tmpout out convergenceMG_history.txt

