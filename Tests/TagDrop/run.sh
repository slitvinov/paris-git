#!/bin/bash
#set -x


/bin/rm -fr out input reference.txt
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
    precision=1e-5
fi

mpirun -np $npstart paris > tmpout 2>&1
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

echo -e "$GREEN" "Check results using Visit/Tecplot."  "$NORMAL"

#if [ -d out ]; then
#    cd out
## assume that curvature.txt was written by paris
#      cat curvature-0000?.txt > curvature.txt
#      cat reference-0000?.txt > reference.txt
#      compare curvature.txt reference.txt $precision
#else
#    RED="\\033[1;31m"
#    NORMAL="\\033[0m"
#    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
#fi


