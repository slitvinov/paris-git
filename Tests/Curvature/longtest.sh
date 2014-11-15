#!/bin/bash
#set -x

d=2
ndepth=3
samplesize=16

here=`pwd`
#runtest.sh 16 $samplesize $d
#cd ../../Devel/Curvature-test
#./compare-inf.sh $ndepth $d
#cd $here

d=3
ndepth=3

here=`pwd`
runtest.sh 16 $samplesize $d || (echo "Failed curvature test";  exit)
cd ../../Devel/Curvature-test

./compare-inf.sh $ndepth $d


