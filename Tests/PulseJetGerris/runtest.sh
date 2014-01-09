#! /bin/bash
#set -x

if [ $# -lt 3 ]; then
    echo "missing arguments"
    echo usage $0 dt nx Implicit:T/F
    exit
fi

dt=`echo $1 | sed s/d/e/g | sed s/D/e/g`   # should always be checked ! 
let nx=$2
imp=$3
maxerr=`awk -v dt=$dt 'BEGIN {print 0.01*dt}' `

npx=2; npy=$npx; npz=$npx

/bin/rm -f flowrates-IMP-$imp*

rm -fr input out stats
sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/g | sed s/IMPTEMP/$imp/g > testinput-$nx.tmp
sed s/NPXTEMP/$npx/g testinput-$nx.tmp |  sed s/NPYTEMP/$npy/g |  sed s/NPZTEMP/$npz/g | sed s/MAXERRTEMP/$maxerr/g > testinput-$nx    
ln -s testinput-$nx input
let np=$npx*$npx*$npx
 mpirun -quiet -np $np paris  > tmpout-$nx
#    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout-$nx





