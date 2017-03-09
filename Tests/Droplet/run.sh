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
	nfilter=1
    else
	echo "usage: $0" 
	echo "or:    $0 MomCons"
	exit 1
    fi
fi

sed s/MOMCONSTEMP/$momconstrue/g inputvof.template | sed s/NFILTERTEMP/$nfilter/g > inputvof 
sed s/MOMCONSTEMP/$plotlabel/g plot.gp.template > plot.gp

if [ "$HAVE_SILO" == 1 ]; then
  echo "we have silo"
  sed s/"output_format = 2"/"output_format = 5"/g inputshort > tmp
  mv inputshort inputshort.bkp
  mv tmp inputshort
else
  echo "we do not have silo"
  sed s/"output_format = 5"/"output_format = 2"/g inputshort > tmp
  mv inputshort inputshort.bkp
  mv tmp inputshort
fi

rm -fR input out
ln -s inputshort input
mpirun -np 8 paris > tmpout
rm -f input
rm inputshort
mv inputshort.bkp inputshort

echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
gnuplot < plot.gp 

precision=2e-2
pariscompare out/droplet-test-vel.txt reference-short-$plotlabel.txt $precision 1 1


GREEN="\\033[1;32m"
NORMAL="\\033[0m"



