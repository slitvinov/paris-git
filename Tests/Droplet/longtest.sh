#!/bin/bash
#set -x

./run_one_test.sh 32 2 0.02 5d-6

gnuplot < plot.gp
mv droplet.png ../Testreport

