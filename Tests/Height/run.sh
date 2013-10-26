#!/bin/bash
#set -x

type=height_test

normup=T
cyldir=3

radius=0.32
precision=1e-20
refinement=1

nx=16
ny=$nx; nz=$ny
npx=2; npy=2; npz=2

if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
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


if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	npy=1; npz=1; npx=1
    else
	echo "unknown option" $1 
	exit 1
    fi
fi

let npstart=$npx*$npy*$npz

sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/g | sed s/NZTEMP/$nz/g | sed s/NPZTEMP/$npz/g  | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g > testinput
sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
ln -s testinput-$dim-$nx-$radius input
sed s/REFINEMENTTEMP/$refinement/g inputvof.template | sed s/TYPETEMP/$type/g  | sed s/CYLDIRTEMP/$cyldir/g | sed s/NORMUPTEMP/$normup/g > inputvof

mpirun -np $npstart paris > tmpout 2>&1
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

if [ -d out ]; then
    cd out
	cat height-0000?.txt >> output1
	cat reference-0000?.txt >> reference.txt
	compare output1 reference.txt $precision
    cd ..
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi

