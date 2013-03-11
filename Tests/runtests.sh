#!/bin/bash
#set -x

for dir in `ls`; do 
    if [ -d $dir ]; then
	cd $dir
	if [ -e 'run.sh' ]; then
	    echo running test in $dir
	    sh run.sh > outtest
# last line in output of test should be PASS or FAIL
	    tail -n 1 outtest
	fi
	cd ..
    fi
done
