#! /bin/bash
#set -x

./run_one_test.sh F Curvature2D T 2 0.32 4e-2 8 16
./run_one_test.sh F Curvature2D T 3 0.32 3e-1 8 16
./run_one_test.sh T Curvature2D F 2 0.2 3.6e-1 4 8
./run_one_test.sh T Curvature2D F 3 0.2 3e-1 4 8
echo "large mono"
./run_one_test.sh T Curvature2D F 3 0.32 3e-1 8 16
echo "convergence test"
runtest.sh 8


#ismono=$1
#type=$2
#normup=$3
#cyldir=$4
#radius=$5
#precision=$6
#refinement=$7
#nx=$8
