#! /bin/bash
#set -x



if [ $# -lt 1 ]; then
    echo "missing arguments"
#    echo usage $0 dimension initialisation-over-refinement-level
    echo usage $0 initialisation-over-refinement-level
    exit
fi

#let d=$1
let d=2
dim="$d"D
init=$1 # !!

if [ $d == 2 ]; then
list="3 4 5"
npz=1
nz=2  # thus cyldir=3
cyldir=3
else
list="3 4 5 6 7"
npz=2
fi
echo $list

/bin/rm -f *.tmp
for level in $list; do
    echo $level
    nx=`awk -v level=$level 'BEGIN {print 2**level}'`
    refinement=`awk -v init=$init 'BEGIN {print 2**init}'`
    if [ $d == 3 ]; then
	nz=$nx
    fi

    for radius in 0.2 0.25 0.32; do 
	./run_one_test.sh F Curvature2D F $cyldir $radius 1e20 $init $nx > /dev/null
	awk -v nx=$nx -v radius=$radius '{print nx * radius , $1, $2, $3 }' mcount.tmp >> allcount.tmp
	cd out
	compare curvature.txt reference.txt 1e20 1 2 >> ../cmpout.tmp
	echo `awk -v nx=$nx -v radius=$radius 'BEGIN {print nx * radius }'`  `compare curvature.txt reference.txt 0.1 1 2 `  >> ../out.tmp
	cd ..
    done
done

# one more time for radius 0.4 at max level

radius=0.4
echo $level

./run_one_test.sh F Curvature2D F $cyldir $radius 1e20 $init $nx > /dev/null
cd out
compare curvature.txt reference.txt 1e20 1 2 >> ../cmpout.tmp
echo `awk -v nx=$nx -v radius=$radius 'BEGIN {print nx * radius }'`  `compare curvature.txt reference.txt 0.1 1 2 `  >> ../out.tmp
cd ..



gnuplot <<EOF
set log x
set log y
set xlabel "Grid points per Radius"
set ylabel "Curvature Error"
plot "out.tmp" u 1:2 t "L2", 2/(x*x)
set term pdf
set out "curvature.pdf"
plot "out.tmp" u 1:2 t "L2", 2/(x*x) 
EOF
