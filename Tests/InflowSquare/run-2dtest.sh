#! /bin/bash
#set -x
export LANG=en_EN

if [ $# -lt 5 ]; then
    echo "missing arguments"
    echo usage $0 nr_of_dt_values dt0 nx0 Implicit:T/F precision
    exit
fi

let ndt=$1
dt0=$2
let nx0=$3
let idt=0
imp=$4
precision=$5

dt=$dt0

/bin/rm -f flowrates-IMP-$imp*

# here=`pwd`
# there=../Disk/out
# if [ -f $there/bitmap-00000.txt ]; then
#     cd $there
#     cp bitmap-*.txt $here
#     cd $here
# else
#     echo "error no bitmap-*.txt"
#     exit
# fi

while [ $idt -lt $ndt ] ; do
    rm -fr input out stats
    let nx=$nx0
    zlength=`awk -v nx=$nx 'BEGIN { print 2./nx}'`
    sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/ | sed s/IMPTEMP/$imp/ | sed s/ZLTEMP/$zlength/  > testinput-$nx-$idt
    ln -s testinput-$nx-$idt input
    if [ `grep DoFront input | awk '{print $3}'` == 'T' ]; then
	let np=17
    else
	let np=16
    fi
    mpirun -quiet -np $np paris > tmpout-$nx-$idt
    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout-$nx-$idt
    awk -v dt=$dt '{ print dt " " $1 " " $2}' < out/flowrate.txt >> flowrates-IMP-$imp.txt
    dt=`awk -v dt=$dt 'BEGIN {print dt/2}'`
    let idt=$idt+1
done
 
if [ -d out ]; then
    cd out
	compare ../reference.txt flowrate.txt $precision
    cd ..
else
    echo "FAIL: directory out not created"
fi



