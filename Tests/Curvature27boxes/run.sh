#!/bin/bash
#set -x

./run_one_test.sh F 1 F 1 1 4e-3 1 36 1 1 1 3
pariscompare out/curvature.txt out/reference.txt 4e-3 1 1
./run_one_test.sh T 1 F 1 1 1e-1 1 36 1 1 1 3
pariscompare out/curvature.txt out/reference.txt 7e-3 1 1

