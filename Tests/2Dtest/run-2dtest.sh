#! /bin/bash
#set -x

if [ $# -lt 6 ]; then
    echo "missing arguments"
    echo usage $0 nr_of_dt_values nr_of_nx_values dt0 nx0 Implicit:T/F precision
    exit
fi

let ndt=$1
let ndx=$2
dt0=$3
let nx0=$4
let idt=0
imp=$5
precision=$6

dt=$dt0

/bin/rm flowrates-IMP-$imp*

while [ $idt -lt $ndt ] ; do
    rm -fr input out stats
    let nx=$nx0
    sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/g | sed s/IMPTEMP/$imp/g > testinput-$nx-$idt
    ln -s testinput-$nx-$idt input
    if [ `grep DoFront input | awk '{print $3}'` == 'T' ]; then
	let np=5
    else
	let np=4
    fi
    mpirun -np $np paris > tmpout-$nx-$idt
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
#! /bin/bash


