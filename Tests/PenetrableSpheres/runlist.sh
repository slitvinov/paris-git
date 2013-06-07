#!/bin/bash
#set -x

if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0 dt0 nx0 
    exit
fi

dt0=$1
let nx0=$2
imp=F
dt=$dt0
permfile=perm-$nx0.txt

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

    Tend=`awk -v n=$nrofcenters 'BEGIN {R=0.0625; phi=exp(-4*3.14157*R*R*R*n/3.); print 400*R**2*phi/(54*log(phi)**2) }'`
    echo $Tend
    rm -fr input out stats
    let nx=$nx0
    sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/g | sed s/IMPTEMP/$imp/g  | sed s/ENDTIMETEMP/$Tend/g > testinput-$nx
    ln -s testinput-$nx input
    let npx=`grep -i NPX input |  awk 'BEGIN {FS = "="}{print $2}' | awk '{print $1}'`
    let npy=`grep -i NPY input |  awk 'BEGIN {FS = "="}{print $2}' | awk '{print $1}'`
    let npz=`grep -i NPZ input |  awk 'BEGIN {FS = "="}{print $2}' | awk '{print $1}'`
    if [ `grep DoFront input | awk '{print $3}'` == 'T' ]; then
	let np=$npx*$npy*$npz+1
    else
	let np=$npx*$npy*$npz
    fi
    mpirun -np $np paris > tmpout-$nx-$nrofcenters
    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout-$nx-$nrofcenters
    awk -v dt=$dt '{ print dt " " $1 " " $2}' < out/flowrate.txt >> flowrates-IMP-$imp.txt
    phi=`cat out/porosity.txt | awk '{print $1}'`
    plot.sh
    awk -v phi=$phi -v nr=$nrofcenters '{radius=0.0625; print nr " " phi " " $2/(radius*radius) " "  phi/(54*log(phi)**2) }' < out/flowrate.txt >> $permfile
done
