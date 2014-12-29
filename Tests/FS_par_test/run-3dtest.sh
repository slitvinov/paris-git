#! /bin/bash
#set -x
export LANG=en_EN

if [ $# -lt 4 ]; then
    echo "missing arguments"
    echo usage $0 nPx nPy nPz nx 
    exit
fi

npx=$1
npy=$2
npz=$3
nx=$4

if [[ $(($nx % $npx)) -ne 0 ]] || [[ $(($nx % $npy)) -ne 0 ]] || [[ $(($nx % $npz)) -ne 0 ]] ; then
    echo "Number of grid points per coordinate direction must be an integer multiple of number of procs in that direction." 
    echo usage $0 nPx nPy nPz nx
    exit 
fi

rm -fr input out stats tmpout-* testinput-* div_type.txt RP_out RK_int_RP*
let nprocs=$1*$2*$3;

sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/ | sed s/NPYTEMP/$npy/ | sed s/NPZTEMP/$npz/> testinput-$npx-$npy-$npz-$nx
ln -s testinput-$npx-$npy-$npz-$nx input

let np=$nprocs

mpirun -np $np paris > tmpout-$npx-$npy-$npz-$nx

if [[ $nprocs == 1 ]] ;  then
	awk '{print $1, $11}' stats > cav_vol-serial
else
	awk '{print $1, $11}' stats > cav_vol-$npx-$npy-$npz-$nx
	pariscompare cav_vol-$npx-$npy-$npz-$nx cav_vol-serial 1e-12
fi
