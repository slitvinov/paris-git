#!/bin/bash
#set -x

for dir in `ls`; do 
    if [ -d $dir ]; then
	if ! [ -a $dir/DONOTRUN ] ; then
	    cd $dir
	    if [ -e 'longtest.sh' ]; then
		echo running long test in $dir
		chmod +x longtest.sh
		./longtest.sh 
	    fi
	    cd ..
	fi
    fi
done
