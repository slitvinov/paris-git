#!/bin/bash
#set -x

if [ $# -lt 7 ]; then
    echo "missing arguments"
    echo usage $0 ismono type cyldir radius initialisation-over-refinement-level nx surf_tension
    exit
fi

ismono=$1
type=$2

cyldir=$3

radius=$4
refinement=$5

nx=$6
stension=$7

if [ $ismono == T ]; then
    setmono=mono
    echo "mono"
else
    setmono=parallel
fi

ny=$nx; nz=$ny
npx=2; npy=1; npz=$npx

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

/bin/rm -fr out input capwave* testinput-* *.visit div_type.txt stats tmpout


if [ $setmono == mono ]; then
    npy=1; npz=1; npx=1
fi

let npstart=$npx*$npy*$npz
ylength=`awk -v nx=$nx 'BEGIN { print 2./nx}'`
sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/g | sed s/NZTEMP/$nz/g | sed s/NPZTEMP/$npz/g  | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g | sed s/LTEMP/$ylength/ | sed s/TENSION/$stension/> testinput
sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
ln -s testinput-$dim-$nx-$radius input
sed s/REFINEMENTTEMP/$refinement/g inputvof.template | sed s/TYPETEMP/$type/g  | sed s/CYLDIRTEMP/$cyldir/g > inputvof

mpirun -np $npstart paris > tmpout


