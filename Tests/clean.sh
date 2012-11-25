#!/bin/bash
#set -x


for dir in `ls`; do 
    cd $dir
    echo "in dir " $dir
    echo "rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.*"
    rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.*
    cd ..
done
