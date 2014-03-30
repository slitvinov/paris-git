#!/bin/bash
#set -x

rm -fR out
mpirun -np 8 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

RED="\\033[1;31m"
GREEN="\\033[1;32m"
NORMAL="\\033[0m"

if [ `grep -c Warning tmpout` -gt 0 ]; then
    tail -3 tmpout
    echo -e "$RED" FAIL 
    echo -e "$NORMAL"
else
    echo -e "$GREEN" PASS 
    echo -e "$NORMAL"
fi




