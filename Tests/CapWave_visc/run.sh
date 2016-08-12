#!/bin/bash
#set -x

tmp=`mktemp -d 2>/dev/null || mktemp -d -t 'tmp'`

rm -fR out
cp inputshort input
mpirun -np 6 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

midpoint=$(cat input | awk -F '=' ' /XLENGTH/ {printf("%15.14f", $2*1.5)}' | tr -d ' ')
echo midpoint is $midpoint
ampini=$(head -n 1 interface.dat | awk -v midpoint=$midpoint '{printf("%15.14f",$2 - midpoint)}')
awk -v midpoint=$midpoint '{printf("%15.14f %15.14f\n",$1, ($2-midpoint)/'$ampini'*0.01)}' interface.dat > simu


