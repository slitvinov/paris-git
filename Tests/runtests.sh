#!/bin/bash
#set -x

cp beginreport.html shorttestreport.html
hash gnuplot 2>/dev/null || { echo >&2 "PARIS testing works better with gnuplot but it is not installed."; }

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
		if [ -e 'report.html' ]; then
		    cat report.html >> ../shorttestreport.html
		fi
	    fi
	    cd ..
	fi
    fi
done
cat endreport.html >> shorttestreport.html
hash Open  2>/dev/null || { exit 0; }
Open shorttestreport.html