#!/bin/bash
#set -x

momconstrue=T
nfilter=1
plotlabel=MomCons


if [ $# -gt 0 ]; then
    if [ $1 == nonMomCons ]; then
	echo "using nonMomCons"
	momconstrue=F
	plotlabel=nonMomCons
	nfilter=1
    else
	echo "usage: $0" 
	echo "or:    $0 nonMomCons"
	exit 1
    fi
fi


sed s/MOMCONSTEMP/$momconstrue/g inputvof.template | sed s/NFILTERTEMP/$nfilter/g > inputvof 
sed s/MOMCONSTEMP/$plotlabel/g plot.gp.template > plot.gp



rm -fR out input
ln -s inputlong input
mpirun -np 8 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
rm -f input

precision=1e-2
pariscompare out/droplet-test-vel.txt reference-$plotlabel.txt $precision 1 0

gnuplot < plot.gp

mv droplet.png ../Testreport
GREEN="\\033[1;32m"
NORMAL="\\033[0m"



