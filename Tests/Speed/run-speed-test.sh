#! /bin/bash
#set -x

if [ $# -lt 1 ]; then
    echo "missing list of np"
    exit
fi

# 14 seconds
make mono64-8

# 21 seconds
make mono64-64

    mpirun -np $np paris > tmpout-$size

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
END { print "size= " size "  cores= " cores "  Z/np= " nrofcells*(time2 - time1)/(cpu2-cpu1)/cores "  Z= " nrofcells*(time2 - time1)/(cpu2-cpu1)} ' < tmpout-$size
done

