#!/bin/bash
#set -x

moftrue=F
nfilter=1

if [ $# -gt 0 ]; then
    if [ $1 == MoF ]; then
	echo "using Mof"
	moftrue=T
	nfilter=0
    fi
fi


sed s/MOFTEMP/$moftrue/g inputvof.template | sed s/NFILTERTEMP/$nfilter/g > inputvof 


rm -fR out
mpirun -np 8 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

gnuplot <<EOF
set xlabel "time"
set ylabel "w(0,0,R)"
set term pdf
set out "tmp.pdf"
plot "out/droplet-test-vel.txt" notitle w l, 0 notitle
set term png
set out "tmp.png"
replot
EOF

precision=2e-2
compare out/droplet-test-vel.txt reference.txt $precision 1 1

GREEN="\\033[1;32m"
NORMAL="\\033[0m"



