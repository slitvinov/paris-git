#! /bin/bash
# set -x



if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0 dimension initialisation-over-refinement-level
    exit
fi

let d=$1
dim="$d"D
init=$2
#let iend=2+3-$d
iend=4
echo iend=$iend

if [ $d == 2 ]; then
list="3 4 5"
else
list="3 4 5"
fi
echo $list
/bin/rm gerris-$dim.txt
for level in $list; do
    echo $level
    for radius in 0.2 0.25 0.32; do 
	gerris$dim -m -DIEND=$iend -DINIT=$init -DLEVEL=$level -DRADIUS=$radius -DD=$d curvature.gfs | awk  -v level=$level -v radius=$radius '{print radius*2**level, $3, $6}' >> gerris-$dim.txt
    done
done
radius=0.4
echo $level
gerris$dim -m  -DIEND=$iend -DINIT=$init -DLEVEL=$level -DRADIUS=$radius -DD=$d curvature.gfs | awk  -v level=$level -v radius=$radius '{print radius*2**level, $3, $6}' >> gerris-$dim.txt

#if [ $d == 2 ]; then
#    gfsview2D curvature--1-3-0.2-2.gfs curvature.gfv &
#fi

gnuplot <<EOF
set log x
set log y
set xlabel "Grid points per Radius"
set ylabel "Curvature Error in $dim"
plot "gerris-$dim.txt" u 1:2 t "L2", "gerris-$dim.txt" u 1:3 t "Linfty"
set term pdf
set out "curvature-$dim.pdf"
plot "gerris-$dim.txt" u 1:2 t "L2", "gerris-$dim.txt" u 1:3 t "Linfty"
EOF