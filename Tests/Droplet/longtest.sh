#!/bin/bash
#set -x

momconstrue=T
nfilter=0
sed s/MOMCONSTEMP/$momconstrue/g inputvof.template | sed s/NFILTERTEMP/$nfilter/g > inputvof 

rm -fR out input
ln -s inputlong input
mpirun -np 8 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
rm -f input

precision=1e-2
pariscompare out/droplet-test-vel.txt reference-MomCons.txt $precision 1 0

gnuplot < plotlong.gp
GREEN="\\033[1;32m"
NORMAL="\\033[0m"



