#!/bin/bash
#set -x



if [ $# -lt 1 ]; then
    echo "missing arguments"
    echo usage $0 endtime
    exit
fi

endtime=$1

rm -fR input out
sed s/ENDTEMP/$endtime/g testinput.template > testinput.tmp
ln -s testinput.tmp input

mpirun -np 8 paris > tmpout

echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

radius=1.5e-3
u0=8.046
pi=3.14157
rho2=1e3
awk -v radius=$radius -v u0=$u0 -v pi=$pi -v rho2=$rho2 '{ print $1, $14/(0.5*u0*u0*rho2*radius*radius*radius*4.0*pi/3.0)}' stats | tail -n +2 > E-k2.tmp
gnuplot < plot.gp 

# use a very lax tolerance for the comparison, since the energy strongly depends on the Poisson solver tolerance
# and we jusst want to see if the code blows up. 
precision=0.1
pariscompare E-k2.tmp reference-E-k2.txt $precision 0 1


