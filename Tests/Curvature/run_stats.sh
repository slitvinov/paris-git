#! /bin/bash
#set -x

if [ $# -lt 4 ]; then
    echo "missing arguments"
    echo usage $0 nx init radius nstats
    exit
fi

nx=$1
init=$2
radius=$3
ny=$nx; nz=$nx
npx=2; npy=$npx; npz=$npx
nstats=$4

let np=$npx*$npy*$npz

ndepth=`head -50  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
dim=3D
/bin/rm -rf errors.tmp

for i in $(seq 1 $nstats); 
do 
echo -n "."
xc=`awk -v nx=$nx -v random=$RANDOM 'BEGIN {print (0.5 + (random / nx)/32767)}'`
yc=`awk -v nx=$nx -v random=$RANDOM 'BEGIN {print (0.5 + (random / nx)/32767)}'`
zc=`awk -v nx=$nx -v random=$RANDOM 'BEGIN {print (0.5 + (random / nx)/32767)}'`

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
    echo `awk -v nx=$nx -v radius=$radius 'BEGIN {print nx * radius }'`  `compare curvature.txt reference.txt 0.1 1 2 `  >> ../errors.tmp
    cd ..
    awk -v nx=$nx  -v radius=$radius '{print nx * radius, $1, $2, $3 }' mcount.tmp >> method_count-$nx-$ndepth.tmp
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
    exit
fi
done
echo " "

cat errors.tmp | awk 'BEGIN {err=0.; a=0.} { if ($3 > err) {err = $3}; a = a +$2}; END {print $1, a/NR,  err, NR}' 
cat errors.tmp | awk 'BEGIN {err=0.; a=0.} { if ($3 > err) {err = $3}; a = a +$2}; END {print $1, a/NR,  err}'  >> paris-$ndepth.tmp

