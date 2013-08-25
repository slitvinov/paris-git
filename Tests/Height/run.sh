#!/bin/bash
#set -x


/bin/rm -fr out input
let npstart=4
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
	cat height-0000?.txt >> output1
	cat reference-0000?.txt >> reference.txt
	compare output1 reference.txt $precision
    cd ..
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi

