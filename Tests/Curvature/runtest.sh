#! /bin/bash
#set -x

if [ $# -lt 4 ]; then
    echo "missing arguments"
    echo usage $0 initialisation-over-refinement-level nstats dim levelmax
    exit 1
fi

init=$1
nstats=$2
d=$3
levelmax=$4

dim=$3'D'
list=`seq 5 $levelmax`
echo $list

ndepth=`head -60  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
/bin/rm -f paris-$dim-$ndepth.tmp method*.tmp

for level in $list; do
    echo $level
    nx=`awk -v level=$level 'BEGIN {print 2**level}'`
    for radius in 0.0625 0.03125 0.1 0.125 0.16 0.2 0.25 0.32; do 
	run_stats.sh $nx $init $radius $nstats $d || exit 1
    done
done

# one more time for radius 0.4 at max level

radius=0.4
echo $level
run_stats.sh $nx $init $radius $nstats $d || exit 1

exit 0 
