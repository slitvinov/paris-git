#!/bin/bash
#set -x

d=3
dim=$d'D'

echo "dimension" $dim
precision=4e-2
./run_one_test.sh F Curvature_test F 0 0.32 4e-2 8 32 0.5 0.5 0.5 3
echo `awk ' /START:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
cd out
pariscompare curvature.txt reference.txt $precision 1 1
cd ..


d=2
dim=$d'D'

echo "dimension" $dim
precision=4e-2
./run_one_test.sh F Curvature2D F 3 0.32 4e-2 8 64 0.5 0.5 0.5 2
echo `awk ' /START:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
cd out
pariscompare curvature.txt reference.txt $precision 1 1
cd ..



#ismono=$1
#type=$2
#normup=$3
#cyldir=$4
#radius=$5

#refinement=$7
#nx=$8
#xc=$9
#yc=$10
#zc=$11
#d=$12
