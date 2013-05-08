#!/bin/bash
#set -x


#echo running test in `pwd`

let np=8
rm -fR out input
ln -s miniinput input
if [ `grep -i DoFront input | awk '{print $3}'` == 'T' ]; then
    let np=$np+1
fi

mpirun -np $np paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

if [ `grep Step tmpout | tail -n 1 |  awk  '{print $2}'` == '3' ]; then 
    echo -e "\033[32;1m PASS\033[0m"
else 
    echo -e "\033[31;1m FAIL\033[0m"
fi




