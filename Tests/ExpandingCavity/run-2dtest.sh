#! /bin/bash
#set -x
export LANG=en_EN

if [ $# -lt 3 ]; then
    echo "missing arguments"
    echo usage $0 dt nx npx
    exit
fi

dt=$1
nx=$2
np_xy=$3

rm -fr input out stats testinput-* *.tmp
zlength=`awk -v nx=$nx 'BEGIN { print 2./nx}'`
let nprocs=$3*$3;
let nprocsfront=$3*$3+1;
sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/ | sed s/ZLTEMP/$zlength/ | sed s/NPXTEMP/$np_xy/ > testinput.tmp
ln -s testinput.tmp input
let np=$nprocs

mpirun -np $np paris > tmpout

awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout
awk -f awk_test.awk stats




