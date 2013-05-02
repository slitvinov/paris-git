#!/bin/bash
#set -x

rm -fR out
mpirun -np 9 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

GREEN="\\033[1;32m"
NORMAL="\\033[0m"

echo -e "$GREEN" "Check results using Visit/Paraview."  "$NORMAL"
