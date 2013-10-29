#!/bin/bash
#set -x

let nx=32
#let nx=16
let nz=$nx
let npx=2
let npz=$npx
radius=0.32
#radius=0.2
refinement=8
dim=3D
precision=4e-2


/bin/rm -fr out input reference.txt

if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	echo "mono"
	let npstart=1
	precision=0.11
	npx=1; npz=1
	radius=0.25
    else
	echo "unknown option" $1 
	exit 1
    fi
fi

let npstart=$npx*$npx*$npz

sed s/NXTEMP/$nx/g testinput.template | sed s/NZTEMP/$nz/g | sed s/NPXTEMP/$npx/g  | sed s/NPZTEMP/$npz/g > testinput
sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
ln -s testinput-$dim-$nx-$radius input
sed s/REFINEMENTTEMP/$refinement/g inputvof.template > inputvof

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



