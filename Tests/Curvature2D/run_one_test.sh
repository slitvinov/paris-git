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
    xc=`awk -F '=' ' /xyzrad\(1, 1\)/ {print $2}' < testinput.template | awk '{print $1}'`
    yc=`awk -F '=' ' /xyzrad\(2, 1\)/ {print $2}' < testinput.template | awk '{print $1}'`
    if [ -f gridgp.template ]; then
	sed s/XC1/$xc/g gridgp.template | sed s/XC2/$yc/g > grid.gp
    fi
fi
if [ $cyldir -gt 3 ]; then
    echo "incorrect cyldir"
    exit
fi

if [ $setmono == mono ]; then
    npy=1; npz=1; npx=1
fi

/bin/rm -fr out input

let npstart=$npx*$npy*$npz

sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/g | sed s/NZTEMP/$nz/g | sed s/NPZTEMP/$npz/g  | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g > testinput
sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
ln -s testinput-$dim-$nx-$radius input
sed s/REFINEMENTTEMP/$refinement/g inputvof.template | sed s/TYPETEMP/$type/g  | sed s/CYLDIRTEMP/$cyldir/g | sed s/NORMUPTEMP/$normup/g > inputvof > inputvof

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

if [ $setmono == mono ] && [ $cyldir == 3 ]; then
	gnuplot <<EOF
call "../grid.gp"
EOF
fi
