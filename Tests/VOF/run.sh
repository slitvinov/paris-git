#!/bin/bash
#set -x

rm -fR out
mpirun -np 8 paris
echo "Check results using Visit/Paraview."

