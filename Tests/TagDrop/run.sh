#!/bin/bash
#set -x


/bin/rm -fr out input stats *root tmpout vol-list.txt  
let npstart=64
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
    precision=1e-5
fi

mpirun -np $npstart paris > tmpout 2>&1
echo `awk ' /START:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`


GREEN="\\033[1;32m"
NORMAL="\\033[0m"

if [ -d out ]; then
      awk '{print $10}' out/element-stats_00000.dat | sort -g | tail -n 4 > vol-list.txt
      pariscompare vol-list.txt ref-vol-list.txt $precision
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi


