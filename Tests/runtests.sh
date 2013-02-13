#!/bin/bash
#set -x

for dir in `ls`; do 
    if [ -d $dir ]; then
	cd $dir
	if [ -e 'run.sh' ]; then
	    echo running test in $dir
	    sh run.sh > outtest
	fi
	cd ..
    fi
done
