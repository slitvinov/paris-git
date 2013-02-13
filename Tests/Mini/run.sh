#!/bin/bash
#set -x


echo running test in `pwd`
rm -fR out input
ln -s miniinput input
mpirun -np 8 paris
