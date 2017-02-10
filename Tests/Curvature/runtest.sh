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

# TUNABLE PARAMETER minlevel (grid is (2**$level)^d as in gerris)

minlevel=5
ndepth=`head -60  ../../surface_tension.f90 |  awk -F '=' ' /NDEPTH/ {print $2}' | tr -d ' '`
# if (2**$minlevel)$npx < 2*ndepth + 3 we will have a problem
# but it will make the script much longer if we test it here. 
# so the next best thing is a crude check

if [ $ndepth -gt 3 ] && [ $minlevel -lt 5 ]; then
    echo "ndepth too large for expected box size nx/npx"
    exit 1
fi

list=`seq $minlevel $levelmax`
echo "list of levels:" $list

/bin/rm -f paris-$dim-$ndepth.tmp method*.tmp

for level in $list; do
    echo "testing level $level"
    nx=`awk -v level=$level 'BEGIN {print 2**level}'`
    if [ $level -eq $minlevel ]; then
	radlist="0.03125 0.04 0.05 0.0625 0.08 0.1 0.125 0.16 0.2 0.25 0.32"
#	radlist="0.03125 0.04 0.05 0.0625 0.08 0.1 0.125 0.16 0.2"
    else
	radlist="0.2 0.25 0.32"
    fi
    for radius in $radlist; do 
	./run_stats.sh $nx $init $radius $nstats $d || exit 1
    done
done

# one more time for radius 0.4 at max level
if [ $minlevel -ne $levelmax ]; then
    radius=0.4
    echo "testing level $level"
    ./run_stats.sh $nx $init $radius $nstats $d || exit 1
else
# I do not remember why - maybe spurious
    echo "not redoing maxlevel"
fi
exit 0 
