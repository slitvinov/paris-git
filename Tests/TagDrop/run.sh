#!/bin/bash
#set -x


/bin/rm -fr out input stats *root tmpout vol-list.txt  
let npstart=8
if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	echo "mono"
	let npstart=1
	ln -s testinput.mono input
	precision=1e-2
    else
	echo "unknown option" $1
	exit 1
    fi
else
    ln -s testinput input
    precision1=1e-5
    precision2=1e-3
fi

mpirun -np $npstart paris > tmpout 2>&1
echo `awk ' /START:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`


GREEN="\\033[1;32m"
NORMAL="\\033[0m"

if [ -d out ]; then
      awk '{print $10}' out/element-stats_00000.dat | sort -g | tail -n 4 > vol-list.txt
      awk '{print $11}' out/element-stats_00000.dat | sort -g | tail -n 4 > sur-list.txt
      pariscompare vol-list.txt ref-vol-list.txt $precision1
      pariscompare sur-list.txt ref-sur-list.txt $precision2
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi


