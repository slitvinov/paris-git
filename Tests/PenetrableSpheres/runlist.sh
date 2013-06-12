#!/bin/bash
#set -x

if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0 dt0 nx0 
    echo "dt0 = 0 for porosity computations"
    exit
fi

dt0=$1
let nx0=$2
imp=F
dt=$dt0
if [ $dt = 0 ]; then
    permfile=phi-$nx0.txt
else
    permfile=perm-$nx0.txt
fi

npx=2
npy=$npx
npz=2

/bin/rm -f $permfile

awk '{print $3}' < centers/centers.txt > tmplist
for nrofcenters in `cat tmplist` ; do
    /bin/rm -f flowrates-IMP-$imp*
    cat < inputsolids.BEGIN > inputsolids
    centers2input.sh centers/centers_$nrofcenters.txt 
    cat < inputcenters.tmp >> inputsolids
    cat >> inputsolids <<EOF
&end
! end of the namelist
EOF
    
    factor=`grep " $nrofcenters " centers/centers.txt | awk  '{ print $7 }' `
    Tend=`awk -v n=$nrofcenters -v f=$factor 'BEGIN {R=0.0625; phi=exp(-4*3.14157*R*R*R*n/3.); print f*400*R**2*phi/(54*log(phi)**2) }'`

    if [ $dt == 0 ]; then
	echo nr of spheres = $nrofcenters
	Tend=0
    else
	echo nr of spheres = $nrofcenters Tend = $Tend factor = $factor
    fi
    rm -fr input out stats
    let nx=$nx0
    sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/g | sed s/IMPTEMP/$imp/g | sed s/ENDTIMETEMP/$Tend/g > testinput-$nx.tmp
    sed s/NPXTEMP/$npx/g testinput-$nx.tmp |  sed s/NPYTEMP/$npy/g |  sed s/NPZTEMP/$npz/g > testinput-$nx
    ln -s testinput-$nx input

    if [ `grep DoFront input | awk '{print $3}'` == 'T' ]; then
	let np=$npx*$npx*$npx+1
    else
	let np=$npx*$npx*$npx
    fi
    mpirun -np $np paris > tmpout-$nx-$nrofcenters
    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout-$nx-$nrofcenters

    phi=`cat out/porosity.txt | awk '{print $1}'`
if [ $dt = 0 ]; then
    awk -v phi=$phi -v nr=$nrofcenters 'BEGIN { printf "%d  %.2f\n", nr, 100*phi}'  >> $permfile
else
    awk -v dt=$dt '{ print dt " " $1 " " $2}' < out/flowrate.txt >> flowrates-IMP-$imp.txt
    plot.sh
    cp tmp.pdf tmp-history-$nrofcenters.pdf
    awk -v phi=$phi -v nr=$nrofcenters '{radius=0.0625; print nr " " phi " " $2/(radius*radius) " "  phi/(54*log(phi)**2) }' < out/flowrate.txt >> $permfile
fi

done
