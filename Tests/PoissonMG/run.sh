#!/bin/bash
#set -x

rm -rf tmpout tmpres

touch tmpres
for it in $(seq 1 1 20); do

rm -rf out
sed 's/MAXITNUM/'$it'/g' inputtemplate > input
mpirun -np 2 paris > tmpout
cpu=$(echo `awk ' /cpu/ { cpu = $3 } END { print cpu } ' < tmpout`)
res=$(echo `awk ' /residual/ { res = $2 } END { print res } ' < tmpout`)
err=$(echo `awk ' /error/ { res = $2 } END { print res } ' < tmpout`)
echo $res $err $cpu >> tmpres
done

echo " "
echo "cpu =" $cpu
awk '{if ('$res' < 0.00000001 && '$err' < 0.1) {print "\033[32;1m PASS\033[0m Residual '$res' Error '$err' "} 
         else {print "\033[31;1m FAIL\033[0m Residual '$res' Error '$err' "}}' tmpres | tail -1

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
plot "tmpres" u 0:1 t 'residual' w lp pt 7, "tmpres" u 0:2 t 'error' w lp pt 7, "tmpres" u 0:3 t 'CPU time' ax x1y2 w p pt 7
pause 3
set term pdf
set out 'convergence.pdf'
replot
exit
EOF

rm -rf input tmpout tmpres out

