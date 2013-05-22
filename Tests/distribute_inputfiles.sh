#! /bin/bash
#set -x

for dir in `ls`
do 
#    echo $dir
    if [ -d $dir ]
    then 
	if [ $dir != 'VOF' ] 
	then 
	    cp VOF/inputvof $dir
	    echo $dir is copied to 
	fi
    fi
 done
