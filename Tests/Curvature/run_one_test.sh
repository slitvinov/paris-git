#!/bin/bash
#set -x


if [ $# -lt 12 ]; then
    echo "missing arguments"
    echo usage $0 ismono type normup cyldir radius precision initialisation-over-refinement-level nx xc yc zc d
    exit
fi

ismono=$1
# type=$2

normup=$3
cyldir=$4

radius=$5
precision=$6
refinement=$7

nx=$8

xc=$9
yc=${10}
zc=${11}
d=${12}

if [ $ismono == T ]; then
    setmono=mono
    echo "mono"
else
    setmono=parallel
fi

ny=$nx; nz=$ny
npx=2; npy=$npx; npz=$npx

if [ $d == 2 ]; then
    type="Curvature2D"
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
	if [ -f grid_template.gp ]; then
	    sed s/XC1/$xc/g grid_template.gp | sed s/XC2/$yc/g  | sed s/RADIUS/$radius/g  > grid.gp
	fi
    fi
    if [ $cyldir -gt 3 ]; then
	echo "incorrect cyldir"
	exit
    fi
else
    npx=2; npy=$npx; npz=$npx
    type="Curvature_test"
    cyldir=0
fi

if [ $setmono == mono ]; then
    npy=1; npz=1; npx=1
fi

/bin/rm -fr out input

ndepth=`head -50  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
dim=$d'D'

let npstart=$npx*$npy*$npz

sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/g | sed s/NZTEMP/$nz/g  | sed s/NPZTEMP/$npz/g  | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g > testinput
sed s/RADIUSTEMP/$radius/g testinput | sed s/XCTEMP/$xc/g  | sed s/YCTEMP/$yc/g  | sed s/ZCTEMP/$zc/g > testinput-$dim-$nx-$radius 
ln -s testinput-$dim-$nx-$radius input

sed s/REFINEMENTTEMP/$refinement/g inputvof.template | sed s/TYPETEMP/$type/g  | sed s/CYLDIRTEMP/$cyldir/g | sed s/NORMUPTEMP/$normup/g > inputvof 

mpirun -np $npstart paris > tmpout 2>&1

success=`tail -1 tmpout | grep -c "Paris exits succesfully after HF"`
if [ $success -ne 1 ]; then
    exit 1
fi
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

if [ $setmono == mono ] && [ $cyldir == 3 ] && [ $nx -le 16 ] && [ d==2 ]; then
	gnuplot <<EOF
call "../grid.gp"
EOF
fi
exit 0
