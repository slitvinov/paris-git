#!/bin/bash
#set -x

let nx=16

if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	echo "mono"
	nx=8
    fi
fi

ny=$nx; nz=$ny
npx=2; npy=2; npz=2

cyldir=3
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

radius=0.32
refinement=4
#type=cylinder_heights
type=Curvature2D

/bin/rm -fr out input reference.txt
let npstart=4
if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	npstart=1
	npy=1
	npz=1
	npx=1
	precision=1e-1
	radius=0.2
    else
	echo "unknown option" $1 
	exit 1
    fi
else
    precision=4e-2
fi
sed s/NXTEMP/$nx/g testinput.template | sed s/NZTEMP/$nz/g | sed s/NPXTEMP/$npx/g  | sed s/NPZTEMP/$npz/g  | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g > testinput
sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
ln -s testinput-$dim-$nx-$radius input
sed s/REFINEMENTTEMP/$refinement/g inputvof.template | sed s/TYPETEMP/$type/g  | sed s/CYLDIRTEMP/$cyldir/g > inputvof

mpirun -np $npstart paris > tmpout 2>&1
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

if [ -d out ]; then
    cd out
      cat curvature-0000?.txt >> curvature.txt
      cat reference-0000?.txt >> reference.txt
      cat bigerror-0000?.txt >> bigerror.txt
      compare curvature.txt reference.txt $precision 1 1
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi

if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	gnuplot <<EOF
call "../grid.gp"
EOF
    fi
fi
