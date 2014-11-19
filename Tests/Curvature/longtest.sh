#!/bin/bash
#set -x

d=2
ndepth=`head -60  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
samplesize=16
levelmax=5

here=`pwd`
runtest.sh 16 $samplesize $d $levelmax || (echo "Failed curvature test";  exit)
cd ../../Devel/Curvature-test
./compare-inf.sh $ndepth $d
cd $here

d=3
here=`pwd`
runtest.sh 16 $samplesize $d  $levelmax || (echo "Failed curvature test";  exit)
cd ../../Devel/Curvature-test
./compare-inf.sh $ndepth $d
cd $here



