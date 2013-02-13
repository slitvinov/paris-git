#!/bin/bash
#set -x

echo running test in `pwd`
rm -fR out
mpirun -np 8 paris
