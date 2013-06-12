#! /bin/bash
#set -x

if [ $# -lt 1 ]; then
    echo "missing arguments"
    echo usage $0 file
    exit
fi


cp $1 permgp.tmp
gnuplot <<EOF &
load "perm.gp"
pause 10
set term pdf
set out 'tmp.pdf'
load "perm.gp"
exit
EOF
