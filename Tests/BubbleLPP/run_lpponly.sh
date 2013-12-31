#!/bin/bash
#set -x


/bin/rm -fr out input inputlpp reference.txt
let npstart=8
if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	   echo "mono"
	   let npstart=1
#      sed 's/out_dir/out_DNS/g' < input_template.mono > input_DNS.mono 
#	   ln -s input_DNS.mono input
#      sed 's/ForT/F/g' < inputlpp_template > inputlpp_DNS 
#	   ln -s inputlpp_DNS inputlpp
#      mpirun -np $npstart paris > tmpout 2>&1
#      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
#      rm  input inputlpp 
      sed 's/out_dir/out_LPP/g' < input_template.mono > input_LPP.mono 
	   ln -s input_LPP.mono input
      sed 's/ForT/T/g' < inputlpp_template > inputlpp_LPP 
	   ln -s inputlpp_LPP inputlpp
      mpirun -np $npstart paris > tmpout 2>&1
      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
    else
	   echo "unknown option" $1
	   exit 1
    fi
else
#      sed 's/out_dir/out_DNS/g' < input_template > input_DNS 
#	   ln -s input_DNS input
#      sed 's/ForT/F/g' < inputlpp_template > inputlpp_DNS 
#	   ln -s inputlpp_DNS inputlpp
#      mpirun -np $npstart paris > tmpout 2>&1
#      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
#      rm  input inputlpp 
      sed 's/out_dir/out_LPP/g' < input_template > input_LPP 
	   ln -s input_LPP input
      sed 's/ForT/T/g' < inputlpp_template > inputlpp_LPP 
	   ln -s inputlpp_LPP inputlpp
      mpirun -np $npstart paris > tmpout 2>&1
      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
fi


gnuplot <<EOF &
set size square
set grid
set nolabel
set title "Velocity of a single bubble" 
set xlabel "time step (dt=5e-6s)"
set ylabel "v (m/s)"
plot 'out_LPP/element-00001.dat' using 1:6 title "LPP"
pause 4
set term pdf
set out 'bubble_vel_lpp.pdf'
plot 'out_LPP/element-00001.dat' using 1:6 title "LPP"
exit
EOF

echo -e "$GREEN" "PASS & Check results in bubble_vel.pdf."  "$NORMAL"


