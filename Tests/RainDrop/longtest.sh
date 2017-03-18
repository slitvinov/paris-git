#!/bin/bash
#set -x

./run_one_test.sh 0.05

gnuplot < plot.gp
mv raindrop-E-k2.png ../Testreport


