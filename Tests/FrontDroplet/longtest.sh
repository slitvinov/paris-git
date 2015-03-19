#!/bin/bash
#set -x

momconstrue=F
plotlabel=NonMomCons
nfilter=1

if [ $# -gt 0 ]; then
    if [ $1 == MomCons ]; then
	echo "using MomCons"
	momconstrue=T
	plotlabel=MomCons
	nfilter=0
    fi
fi
sed s/MOMCONSTEMP/$momconstrue/g inputvof.template | sed s/NFILTERTEMP/$nfilter/g > inputvof 
sed s/nonMOMCONS/$plotlabel/g plot.gp.template > plot.gp

rm -fR out input
ln -s inputlong input
mpirun -np 9 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
rm -f input

gnuplot <<EOF
set xlabel "time"
set ylabel "w(0,0,R)"
set term png
set out "tmp.png"
plot "out/droplet-test-vel.txt" notitle w l, 0 notitle
EOF

#precision=1e-2
#pariscompare out/droplet-test-vel.txt referencelong.txt $precision 1 0

GREEN="\\033[1;32m"
NORMAL="\\033[0m"



