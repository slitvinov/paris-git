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

mv input input.bkp
cp inputlong input
rm -fR out
mpirun -np 8 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
mv input.bkp input

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

#precision=1e-2
#compare out/droplet-test-vel.txt referencelong.txt $precision 1 0

GREEN="\\033[1;32m"
NORMAL="\\033[0m"



