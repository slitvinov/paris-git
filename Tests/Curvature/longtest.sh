#!/bin/bash
#set -x

d=2
ndepth=`head -60  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
samplesize=16
levelmax=6
do2D=1


if [ $do2D == 1 ] 
then
    echo "Launching "$d"D test"
    here=`pwd`
    sh runtest.sh 16 $samplesize $d $levelmax || { 
	echo "Failed curvature test" 
	exit 
    }
    cd ../../Devel/Curvature-test
    ./compare-only-inf.sh $ndepth $d
    cd $here
fi
dim=$d'D'
if [ -d ../Testreport ] ; then
    mv ../../Devel/Curvature-test/curvature$dim.png ../Testreport # unused for now.
fi

d=3
echo "Launching "$d"D test"
here=`pwd`
sh runtest.sh 16 $samplesize $d $levelmax || { 
echo "Failed curvature test" 
exit 
}
cd ../../Devel/Curvature-test
./compare-only-inf.sh $ndepth $d
cd $here

dim=$d'D'
if [ -d ../Testreport ] ; then
    mv ../../Devel/Curvature-test/curvature$dim.png ../Testreport
fi


sed s/SAMPLESIZE/$samplesize/g report-template.html | sed s/NDEPTH/$ndepth/g  > report.html

