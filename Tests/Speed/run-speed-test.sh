#! /bin/bash
#set -x

if [ $# -lt 1 ]; then
    echo "missing list of np"
    echo "typically $0 1 2 4 8 16 to test up to 4096 procs"
    exit
fi

# input size of your machine here
let maxcores=64
let base=64
# mono64-1  to mono64-512  run in less than 30 seconds on babbage cluster
for rootsize in $* ; do
    rm -fr input out
    let nx=$rootsize*$base
    let size=$rootsize*$rootsize*$rootsize
    sed s/NXTEMP/$nx/g input-mono64-template | sed s/NPXTEMP/$rootsize/g | sed s/NPROCESSES/$size/g  > input-mono64-$size
    ln -s input-mono64-$size input
    if [ `grep -i DoFront input | awk '{print $3}'` == 'T' ]; then
	let np=$size+1
    else
	let np=$size
    fi
    echo "number of procs, np variable " $np
    mpirun -quiet -np $np paris > tmpout-$size

    let Cores=$size
    let nrofcells=$nx*$nx*$nx
    if [ $size -gt $maxcores ]; then
	Cores=$maxcores
    fi
    awk -v nrofcells=$nrofcells -v cores=$Cores -v size=$size ' /Step:/ { 
    if ($2 == 1) {
	time1 = $2;
	cpu1 = $8;
    }
    else if   ($2 == 11) {
	time2 = $2;
	cpu2 = $8;
    }
} 
END { print "\033[1;32m " "size= " size "  cores= " cores "  Z/np= " nrofcells*(time2 - time1)/(cpu2-cpu1)/cores "  Z= " nrofcells*(time2 - time1)/(cpu2-cpu1) " \033[0m" } ' < tmpout-$size
done

