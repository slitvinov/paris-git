#!/bin/bash
#set -x

#!/bin/bash
#set -x

let nx=32
let nz=$nx
let npx=2
let npz=$npx
radius=0.32
refinement=4
dim=3D

/bin/rm -fr out input reference.txt
let npstart=$npx*$npx*$npz
if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	echo "mono"
	let npstart=1
	ln -s testinput.mono input
	precision=5e-2
    else
	echo "unknown option" $1 
	exit 1
    fi
else
    sed s/NXTEMP/$nx/g testinput.template | sed s/NZTEMP/$nz/g | sed s/NPXTEMP/$npx/g  | sed s/NPZTEMP/$npz/g > testinput
    sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
    ln -s testinput-$dim-$nx-$radius input
    sed s/REFINEMENTTEMP/$refinement/g inputvof.template > inputvof
    precision=4e-2
fi

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



