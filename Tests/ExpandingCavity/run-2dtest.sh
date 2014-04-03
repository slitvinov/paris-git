#! /bin/bash
#set -x
export LANG=en_EN

if [ $# -lt 3 ]; then
    echo "missing arguments"
    echo usage $0 dt nx Implicit:T/F 
    exit
fi

dt=$1
nx=$2
np_xy=$3
imp=$4

rm -fr input out stats
zlength=`awk -v nx=$nx 'BEGIN { print 2./nx}'`
let nprocs=$3*$3;
let nprocsfront=$3*$3+1;
sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/ | sed s/IMPTEMP/$imp/ | sed s/ZLTEMP/$zlength/ | sed s/NPXTEMP/$np_xy/ > testinput.tmp
ln -s testinput.tmp input
#   if [ `grep DoFront input | awk '{print $3}'` == 'T' ]; then
#	let np=5
#   else
let np=$nprocs
#  fi
mpirun -np $np paris > tmpout

RED="\\033[1;31m"
GREEN="\\033[1;32m"
NORMAL="\\033[0m"

if [ `grep -c ERROR tmpout ` -gt 0 ]; then
    echo -e "$RED" FAIL "$NORMAL"
else
    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout
    echo -e "$GREEN" PASS "$NORMAL"
exit # exit to avoid graphic output on comput. servers
    ./gnmovie.sh
    ./volumeplot.sh
fi




