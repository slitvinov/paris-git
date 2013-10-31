#! /bin/bash
#set -x

if [ $# -lt 6 ]; then
    echo "missing arguments"
    echo usage $0 nx init radius xc yc zc
    exit
fi

nx=$1
init=$2
radius=$3
xc=$4
yc=$5
zc=$6

ny=$nx; nz=$nx
npx=2; npy=$npx; npz=$npx


let np=$npx*$npy*$npz

ndepth=`head -50  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
dim=3D

sed s/NXTEMP/$nx/g testinput.template | sed s/NZTEMP/$nz/g | sed s/NPXTEMP/$npx/g  | sed s/NPZTEMP/$npz/g > testinput
sed s/RADIUSTEMP/$radius/g testinput | sed s/XCTEMP/$xc/g  | sed s/YCTEMP/$yc/g  | sed s/ZCTEMP/$zc/g > testinput-$dim-$nx-$radius 
sed s/REFINEMENTTEMP/$init/g inputvof.template > inputvof
rm -fr input out stats 
ln -s testinput-$dim-$nx-$radius input
mpirun -np $np paris > tmpout
if [ -d out ]; then
    cd out
    cat curvature-0000?.txt >> curvature.txt
    cat reference-0000?.txt >> reference.txt
    echo `awk -v nx=$nx -v radius=$radius 'BEGIN {print nx * radius }'`  `compare curvature.txt reference.txt 0.1 1 2 `  >> ../paris-$nx-$ndepth.tmp
    cd ..
    awk -v nx=$nx  -v radius=$radius '{print nx * radius, $1, $2, $3 }' mcount.tmp >> method_count-$nx-$ndepth.tmp
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi
