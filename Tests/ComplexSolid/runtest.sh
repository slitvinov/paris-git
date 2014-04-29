#! /bin/bash
#set -x

if [ $# -lt 5 ]; then
    echo "missing arguments"
    echo usage $0 nr_of_dt_values dt0 nx Implicit:T/F precision 
    exit
fi

let ndt=$1
dt0=`echo $2 | sed s/d/e/g | sed s/D/e/g`   # should always be checked ! 
let nx=$3
let idt=0
imp=$4
precision=$5
maxerr=`awk -v dt0=$dt0 'BEGIN {print 0.001*dt0}'`

dt=$dt0

npx=2; npy=$npx; npz=$npx

# /bin/rm -f flowrates-IMP-$imp*

# #prepare rock

if [ -s bitmap-0000.txt ]; then
    rm bitmap*.txt
fi

echo "preparing Rock$nx"
here=`pwd`
cd Rock$nx/
rockread $nx $npx < rock_$nx.txt > rockread_out.tmp
mv bitmap*.txt $here
cd $here


while [ $idt -lt $ndt ] ; do
    rm -fr input out stats
    sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/g | sed s/IMPTEMP/$imp/g > testinput-$nx-$idt.tmp
    sed s/NPXTEMP/$npx/g testinput-$nx-$idt.tmp |  sed s/NPYTEMP/$npy/g |  sed s/NPZTEMP/$npz/g | sed s/MAXERRTEMP/$maxerr/g > testinput-$nx-$idt
    ln -s testinput-$nx-$idt input
    let np=$npx*$npx*$npx
    mpirun -quiet -np $np paris > tmpout-$nx-$idt
    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout-$nx-$idt
    dt=`awk -v dt=$dt 'BEGIN {print dt/2}'`
    let idt=$idt+1
done

#phi=`cat out/porosity.txt | awk '{print $1}'`
#echo phi = $phi
awk '{ print $1,$3 } ' < stats > volume.tmp
awk '{ print $1,$12 } ' < stats > centerofmass.tmp
tail volume.tmp > compare.tmp
cat centerofmass.tmp >> compare.tmp
tail volumeref.txt > compareref.tmp
cat centerofmassref.txt >> compareref.tmp

precision=1e-4
compare  compareref.tmp compare.tmp $precision 0 0


