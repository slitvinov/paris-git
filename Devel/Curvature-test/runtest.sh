#! /bin/bash
#set -x

if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0 dimension initialisation-over-refinement-level
    exit
fi

let d=$1
dim="$d"D
init=$2

if [ $d == 2 ]; then
list="3 4 5 6 7"
npz=1
nz=2
else
list="3 4 5 6 7"
npz=2
fi
echo $list

npx=2
let np=$npx*$npx*$npz


/bin/rm -f *.tmp
for level in $list; do
    echo $level
    nx=`awk -v level=$level 'BEGIN {print 2**level}'`
    refinement=`awk -v init=$init 'BEGIN {print 2**init}'`
    if [ $d == 3 ]; then
	nz=$nx
    fi
    sed s/NXTEMP/$nx/g testinput.template | sed s/NZTEMP/$nz/g | sed s/NPXTEMP/$npx/g  | sed s/NPZTEMP/$npz/g > testinput
    for radius in 0.2 0.25 0.32; do 
	sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
	sed s/REFINEMENTTEMP/$refinement/g inputvof.template > inputvof
	rm -fr input out stats 
	ln -s testinput-$dim-$nx-$radius input
	mpirun -np $np paris > tmpout
	if [ -d out ]; then
	    cd out
	    cat curvature-0000?.txt >> curvature.txt
	    cat reference-0000?.txt >> reference.txt
	    compare curvature.txt reference.txt 1e20 1 2 >> ../cmpout.tmp
    echo `awk -v nx=$nx -v radius=$radius 'BEGIN {print nx * radius }'`  `compare curvature.txt reference.txt 0.1 1 2 `  >> ../paris-$dim.tmp
	    cd ..
	else
	    RED="\\033[1;31m"
	    NORMAL="\\033[0m"
	    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
	fi
    done
done

# one more time for radius 0.4 at max level

radius=0.4
echo $level
sed s/NXTEMP/$nx/g testinput.template | sed s/NZTEMP/$nz/g | sed s/NPXTEMP/$npx/g  | sed s/NPZTEMP/$npz/g > testinput
sed s/RADIUSTEMP/$radius/g testinput > testinput-$dim-$nx-$radius 
rm -fr input out stats 
ln -s  testinput-$dim-$nx-$radius  input
mpirun -np $np paris > tmpout
if [ -d out ]; then
    cd out
    cat curvature-0000?.txt >> curvature.txt
    cat reference-0000?.txt >> reference.txt
    compare curvature.txt reference.txt 1e20 1 2 >> ../cmpout.tmp
    echo `awk -v nx=$nx -v radius=$radius 'BEGIN {print nx * radius }'`  `compare curvature.txt reference.txt 0.1 1 2 `  >> ../paris-$dim.tmp
    cd ..
else
    RED="\\033[1;31m"
    NORMAL="\\033[0m"
    echo -e "$RED" "FAIL: directory 'out' not found."  "$NORMAL"
fi

