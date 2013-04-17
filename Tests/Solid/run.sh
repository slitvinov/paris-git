#!/bin/bash
set -x

rm -fr input out stats tmpout

if [ $# -gt 0 ] ; then
    if [ $1 == "exp" ] ; then
	ln -s testinput.explicit input
	mpirun -np 9 paris > tmpout
    elif [ $1 == "mono" ] ; then
	ln -s testinput.mono input
	paris > tmpout
    else
	echo "$0: invalid option " $1
	exit 1
    fi
else
    ln -s testinput.solid input
    mpirun -np 9 paris > tmpout
fi 


NORMAL="\\033[0;39m"
GREEN="\\033[1;32m"

echo -e "$GREEN" "Check results using Visit/Paraview." "$NORMAL"
