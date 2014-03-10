#! /bin/bash
#set -x
export LANG=en_EN

if [ $# -lt 5 ]; then
    echo "missing arguments"
    echo usage $0 nr_of_dt_values dt0 nx0 Implicit:T/F precision
    exit
fi

let ndt=$1
dt0=$2
let nx0=$3
let idt=0
imp=$4
precision=$5

dt=$dt0

/bin/rm -f flowrates-IMP-$imp*

while [ $idt -lt $ndt ] ; do
    rm -fr input out stats
    let nx=$nx0
    zlength=`awk -v nx=$nx 'BEGIN { print 2./nx}'`
    sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/ | sed s/IMPTEMP/$imp/ | sed s/ZLTEMP/$zlength/  > testinput-$nx-$idt
    ln -s testinput-$nx-$idt input
    if [ `grep DoFront input | awk '{print $3}'` == 'T' ]; then
	let np=5
    else
	let np=4
    fi
    mpirun -np $np paris > tmpout-$nx-$idt
    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout-$nx-$idt
    awk -v dt=$dt '{ print dt " " $1 " " $2}' < out/flowrate.txt >> flowrates-IMP-$imp.txt
    dt=`awk -v dt=$dt 'BEGIN {print dt/2}'`
    let idt=$idt+1
done

end=`grep -i EndTime input |  awk 'BEGIN {FS = "="}{print $2}' | awk '{print $1}'`
phi=`cat out/porosity.txt | awk '{print $1}'`
# echo phi = $phi

awk '{print $1 " " $3}' < stats > deriv.tmp
parisdeconv deriv.tmp > toplot.txt

gnuplot <<EOF > tmp 2>&1  &
set term pdf
set out "tmp.pdf"
f(x) = a*x + b
FIT_LIMIT = 1e-6
fit [2*$end/3:$end] f(x) "toplot.txt" via a, b
plot "toplot.txt", f(x)
print "permeability = ",-1/a
# pause 10
EOF

grep permeability tmp | awk -v phi=$phi '{print $3*phi}' > perm.tmp

awk '{print $2}' < out/flowrate.txt >> perm.tmp
 
if [ -d out ]; then
    cd out
	compare ../reference.txt flowrate.txt $precision
    cd ..
else
    echo "FAIL: directory out not created"
fi
#! /bin/bash


