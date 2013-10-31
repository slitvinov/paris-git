#!/bin/bash
#set -x

precision=4e-2

run_one_test.sh 32 8 0.32 0.5 0.5 0.5

cd out
compare curvature.txt reference.txt $precision 1 1
cd ..




