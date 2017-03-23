#!/bin/bash
#set -x

hash gnuplot 2>/dev/null || { echo >&2 "PARIS testing works better with gnuplot but it is not installed."; }

for dir in `ls`; do 
    if [ -d $dir ]; then
	if ! [ -a $dir/DONOTRUN ] ; then
	    cd $dir
	    if [ -e 'run.sh' ]; then
		echo running test in $dir
		chmod +x run.sh
		./run.sh > outtest 2>&1
# last line in output of test should be PASS or FAIL
		tail -n 2 outtest
	    fi
	    cd ..
	fi
    fi
done
