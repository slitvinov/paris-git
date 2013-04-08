#!/bin/bash
#set -x

echo running test in `pwd`

rm -fR out input
ln -s testinput.solid input

mpirun -np 9 paris

NORMAL="\\033[0;39m"
GREEN="\\033[1;32m"

echo -e "$GREEN" "Check results using Visit/Paraview." "$NORMAL"
