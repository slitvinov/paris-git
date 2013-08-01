#!/bin/bash
#set -x


/bin/rm -fr out input reference.txt
let npstart=8
if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	echo "mono"
	let npstart=1
	ln -s testinput.mono input
	precision=1e-20
    else
	echo "unknown option" $1
	exit 1
    fi
else
    ln -s testinput input
    precision=1e-20
fi

mpirun -np $npstart paris > tmpout 2>&1
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

if [ -d out ]; then
    cd out
# assume that curvature.txt was written by paris
	cat curvature.txt >> ../output
	echo "1 0.22" >> ../reference.txt 
    cd ..
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi

compare output reference.txt $precision

