#!/bin/bash
#set -x

rm -fR out
mpirun -np 9 paris > tmpout

GREEN="\\033[1;32m"
NORMAL="\\033[0m"

echo -e "$GREEN" "Check results using Visit/Paraview."  "$NORMAL"
