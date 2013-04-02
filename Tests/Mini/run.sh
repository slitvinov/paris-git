#!/bin/bash
#set -x


echo running test in `pwd`

let np=8
rm -fR out input
ln -s miniinput input
if [ `grep TWOPHASE input | awk '{print $3}'` == 'T' ]; then
    let np=$np+1
fi

mpirun -np $np paris > tmpout

if [ `grep Step tmpout | tail -n 1 |  awk  '{print $2}'` == '3' ]; then 
    echo "\033[32;1m PASS\033[0m"
else 
    echo "\033[31;1m FAIL\033[0m"
fi
rm tmpout input



