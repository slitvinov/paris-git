#!/bin/bash
#set -x

./run-2dtest.sh 1.0e-4 32 2
awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout
awk -f awk_test.awk stats
