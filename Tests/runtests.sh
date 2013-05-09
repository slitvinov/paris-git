#!/bin/bash
#set -x

for dir in `ls`; do 
    if [ -d $dir ]; then
	if ! [ -a $dir/DONOTRUN ] ; then
	    cd $dir
	    if [ -e 'run.sh' ]; then
		echo running test in $dir
		chmod +x run.sh
		./run.sh > outtest
# last line in output of test should be PASS or FAIL
		tail -n 2 outtest
	    fi
	    cd ..
	fi
    fi
done
