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


rm -fR input out
ln -s inputshort input
mpirun -np 8 paris > tmpout
rm -f input
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

if [ ${CGFONTGETGLYPH_PARIS_PROBLEM:=0} == 0 ]; then
gnuplot <<EOF
set xlabel "time"
set ylabel "w(0,0,R)"
set term png
set out "tmp.png"
load "plot.gp"
EOF
else
gnuplot <<EOF
set xlabel "time"
set ylabel "w(0,0,R)"
set term png
set out "tmp.png"
load "plot.gp"
EOF
fi

precision=2e-2
pariscompare out/droplet-test-vel.txt reference.txt $precision 1 1


GREEN="\\033[1;32m"
NORMAL="\\033[0m"



