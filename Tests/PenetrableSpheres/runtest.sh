#! /bin/bash
#set -x

if [ $# -lt 7 ]; then
    echo "missing arguments"
    echo usage $0 nr_of_dt_values dt0 nx0 Implicit:T/F precision nrofcenters Tend 
    echo "if Tend=0 just compute rock porosity, Tend=-1 automatic computation of Tend"
    echo "if Tend=0 npx=1 otherwise npx=2"
    exit
fi

let ndt=$1
dt0=`echo $2 | sed s/d/e/g | sed s/D/e/g`   # should always be checked ! 
let nx0=$3
let idt=0
imp=$4
precision=$5
nrofcenters=$6
maxerr=`awk -v dt0=$dt0 'BEGIN {print 0.1*dt0}'`

dt=$dt0
if [ $7 == "-1" ]; then
    Tend=`awk -v n=$nrofcenters 'BEGIN {R=0.0625; phi=exp(-4*3.14157*R*R*R*n/3.); print 400*R**2*phi/(54*log(phi)**2) }'`
    echo "using default Tend = " $Tend
else
    Tend=$7
fi

if [ $Tend == 0 ] ; then
    Tend=$dt
    npx=1
else
    npx=2
fi
npy=$npx
npz=1


if ! [ -f centers/centers_$nrofcenters.txt ]; then
  echo -e "\033[31;1m $0: error:  no file centers_$nrofcenters.txt \033[0m"
  exit 1
fi

/bin/rm -f flowrates-IMP-$imp*
cat < inputsolids.BEGIN > inputsolids
centers2input.sh centers/centers_$nrofcenters.txt 
cat < inputcenters.tmp >> inputsolids
cat >> inputsolids <<EOF
&end
! end of the namelist<
EOF

while [ $idt -lt $ndt ] ; do
    rm -fr input out stats
    let nx=$nx0
    sed s/NXTEMP/$nx/g testinput.template | sed s/DTTEMP/$dt/g | sed s/IMPTEMP/$imp/g | sed s/ENDTIMETEMP/$Tend/g > testinput-$nx-$idt.tmp
    sed s/NPXTEMP/$npx/g testinput-$nx-$idt.tmp |  sed s/NPYTEMP/$npy/g |  sed s/NPZTEMP/$npz/g | sed s/MAXERRTEMP/$maxerr/g > testinput-$nx-$idt
    ln -s testinput-$nx-$idt input
    if [ `grep DoFront input | awk '{print $3}'` == 'T' ]; then
	let np=$npx*$npx*$npx+1
    else
	let np=$npx*$npx*$npx
    fi
    mpirun -np $np -quiet paris > tmpout-$nx-$idt
    awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout-$nx-$idt
    awk -v dt=$dt '{ print dt " " $1 " " $2}' < out/flowrate.txt >> flowrates-IMP-$imp.txt
    dt=`awk -v dt=$dt 'BEGIN {print dt/2}'`
    let idt=$idt+1
done

end=`grep -i EndTime input |  awk 'BEGIN {FS = "="}{print $2}' | awk '{print $1}'`
phi=`cat out/porosity.txt | awk '{print $1}'`

if [ $7 == 0 ] ; then 
    echo phi = $phi
    exit
fi

awk '{print $1 " " $3}' < stats > deriv.tmp
parisdeconv deriv.tmp > toplot.tmp

gnuplot <<EOF > tmp 2>&1
f(x) = a*x + b
FIT_LIMIT = 1e-6
fit [2*$end/3:$end] f(x) "toplot.tmp" via a, b
plot "toplot.tmp", f(x)
pause 10
print "permeability = ",-1/a
EOF

grep permeability tmp | awk -v phi=$phi '{print "decaytime = " $3*phi}' > perm.txt

awk '{radius=0.0625; print "K/R^2 = " $2/(radius*radius)}' < out/flowrate.txt >> perm.txt

if [ -d out ]; then
    cd out
	pariscompare ../reference.txt flowrate.txt $precision
    cd ..
else
    echo -e "\033[31;1m FAIL: directory out not created\033[0m"
fi
#! /bin/bash


