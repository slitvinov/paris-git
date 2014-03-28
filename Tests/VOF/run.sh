#!/bin/bash
#set -x

rm -fR out stats
mpirun -np 8 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

awk '{ print $1,$3 } ' < stats > volume.tmp
awk '{ print $1,$12 } ' < stats > centerofmass.tmp
tail volume.tmp > compare.tmp
cat centerofmass.tmp >> compare.tmp
tail volumeref.txt > compareref.tmp
cat centerofmassref.txt >> compareref.tmp


precision=1e-4
compare  compareref.tmp compare.tmp $precision 0 0

GREEN="\\033[1;32m"
NORMAL="\\033[0m"

