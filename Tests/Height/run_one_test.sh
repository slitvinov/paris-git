#!/bin/bash
#set -x

if [ $# -lt 8 ]; then
    echo "missing arguments"
    echo usage $0 ismono type normup cyldir radius precision initialisation-over-refinement-level nx 
    exit
fi

ismono=$1
type=$2

normup=$3
cyldir=$4

radius=$5
precision=$6
refinement=$7

nx=$8

if [ $ismono == T ]; then
    setmono=mono
    echo "mono"
else
    setmono=parallel
fi

ny=$nx; nz=$ny
npx=2; npy=$npx; npz=$npx

if [ $# -gt 0 ]; then
    if [ $setmono == mono ]; then
	echo "mono"
    fi
fi

if [ $cyldir == 1 ];  then
    nx=2
    npx=1
fi
if [ $cyldir == 2 ];  then
    ny=2
    npy=1
fi
if [ $cyldir == 3 ];  then
    nz=2
    npz=1
fi
if [ $cyldir -gt 3 ]; then
    echo "incorrect cyldir"
    exit
fi

/bin/rm -fr out input


if [ $setmono == mono ]; then
    npy=1; npz=1; npx=1
fi

let npstart=$npx*$npy*$npz

sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/g | sed s/NZTEMP/$nz/g | sed s/NPZTEMP/$npz/g  | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g > testinput
sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
ln -s testinput-$dim-$nx-$radius input
sed s/REFINEMENTTEMP/$refinement/g inputvof.template | sed s/TYPETEMP/$type/g  | sed s/CYLDIRTEMP/$cyldir/g | sed s/NORMUPTEMP/$normup/g > inputvof

mpirun -np $npstart paris > tmpout 2>&1
echo `awk ' /START:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

if [ -d out ]; then
    cd out
	cat heighta-0000?.txt >> output1a
	cat heightb-0000?.txt >> output1b
	cat reference-0000?.txt >> reference.txt
	compare output1a reference.txt $precision
	compare output1b reference.txt $precision
    cd ..
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi

