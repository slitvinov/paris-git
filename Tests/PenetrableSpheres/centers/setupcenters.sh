#!/bin/bash
#set -x

if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0 dt0 nx0 
    exit
fi

dt0=$1
let nx0=$2
imp=F
dt=$dt0
file=centers-formatted.txt

/bin/rm -f $file

awk '{print $3}' < centers.txt > tmplist
for nrofcenters in `cat tmplist` ; do
    factor=`grep " $nrofcenters " centers.txt | awk  '{ print $7 }' `
    phipred=`awk -v n=$nrofcenters 'BEGIN {R=0.0625; print 100*exp(-4*3.14157*R*R*R*n/3.)}'`
    phigerris=`grep " $nrofcenters " centers.txt | awk '{print $5}' | tr , .`
    echo $nrofcenters $phipred $phigerris $factor
done
exit
