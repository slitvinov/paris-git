#!/bin/bash
#set -x

rm -fR *~ *.tmp tmp* Testreport*

for dir in `ls`; do 
    if [ -d $dir ]; then
	cd $dir
	echo "Cleaning directory " `pwd`
#	echo "rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.*"
	rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.* statsbub out/*.txt testinput-* *.log aggregate_timer_stats multi*.root error-rank*
	cd ..
    fi
done
