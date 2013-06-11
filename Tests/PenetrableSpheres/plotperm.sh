#! /bin/bash
#set -x

if [ $# -lt 1 ]; then
    echo "missing arguments"
    echo usage $0 nx0 
    exit
fi


cp perm-$1-all.txt permgp.tmp
gnuplot <<EOF &
load "perm.gp"
pause mouse any 
set term pdf
set out 'tmp.pdf'
load "perm.gp"
 exit
EOF
