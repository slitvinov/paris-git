#!/bin/bash
set -x

./run_one_test.sh F Curvature_test F 0 0.25 4e-2 8 32 0.7 0.7 0.5 3
pariscompare out/curvature.txt out/reference.txt 4e-2 1 1
./run_one_test.sh T Curvature_test F 0 0.25 4e-2 8 32 0.7 0.7 0.5 3
pariscompare out/curvature.txt out/reference.txt 4e-2 1 1
./run_one_test.sh F Curvature_test F 0 0.25 4e-2 8 32 0.62 0.62 0.62 3
pariscompare out/curvature.txt out/reference.txt 4e-2 1 1
./run_one_test.sh T Curvature_test F 0 0.25 4e-2 8 32 0.62 0.62 0.62 3
pariscompare out/curvature.txt out/reference.txt 4e-2 1 1
