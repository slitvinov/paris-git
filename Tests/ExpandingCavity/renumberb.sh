#! /bin/bash
#set -x

if [ $# -lt 4 ] 
then
 echo "usage   : " $0 first increment last add [ basename default: goutte ]
 echo "example : " $0 1000  1000      9000  0    nappe
 exit
fi

basename=$5
if [ $# -lt 5 ] 
then
  basename='bubble'
fi

echo 'basename =' $basename

let k=$1
while [ "$k" -le $3 ]
do
    mv $basename.$k.gif $basename.$4$k.gif
    let k+=$2
done
