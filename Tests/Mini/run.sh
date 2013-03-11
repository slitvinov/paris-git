#!/bin/bash
#set -x


echo running test in `pwd`
rm -fR out input
ln -s miniinput input
mpirun -np 8 paris > tmpout
if [ `grep Step tmpout | tail -n 1 |  awk  '{print $2}'` == '100' ]; then 
    echo PASS
else 
    echo FAIL
fi

