#! /bin/bash
#set -x

if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0 initialisation-over-refinement-level nstats
    exit
fi

init=$1
nstats=$2

list="3 4 5 6"
echo $list

ndepth=`head -50  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
/bin/rm -f paris-$ndepth.tmp

for level in $list; do
    echo $level
    nx=`awk -v level=$level 'BEGIN {print 2**level}'`
    for radius in 0.2 0.25 0.32; do 
	run_stats.sh $nx $init $radius $nstats
    done
done

# one more time for radius 0.4 at max level

radius=0.4
echo $level
run_stats.sh $nx $init $radius $nstats
