#!/bin/bash
#set -x


/bin/rm -fr out input inputlpp tmpout*
let npstart=8
if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	   echo "mono"
	   let npstart=1
      sed 's/out_dir/out_LPP/g' < input_template.mono > input_LPP.temp 
      sed 's/dtTBA/2.0e-7/g' < input_LPP.temp > input_LPP.temp1
      sed 's/sigmaTBA/0.0001/g' < input_LPP.temp1 > input_LPP.mono
      rm input_LPP.temp*
	   ln -s input_LPP.mono input
      sed 's/ForT/T/g' < inputlpp_template > inputlpp_LPP 
	   ln -s inputlpp_LPP inputlpp
      mpirun -np $npstart paris > tmpout_LPP 2>&1
      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout_LPP`
      rm  input inputlpp 
      sed 's/out_dir/out_DNS/g' < input_template.mono > input_DNS.temp 
      sed 's/dtTBA/5.0e-8/g' < input_DNS.temp > input_DNS.temp1
      sed 's/sigmaTBA/0.0001/g' < input_DNS.temp1 > input_DNS.mono
      rm input_DNS.temp*
	   ln -s input_DNS.mono input
      sed 's/ForT/F/g' < inputlpp_template > inputlpp_DNS 
	   ln -s inputlpp_DNS inputlpp
      mpirun -np $npstart paris > tmpout_DNS 2>&1
      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout_DNS`
    else
	   echo "unknown option" $1
	   exit 1
    fi
else
      sed 's/out_dir/out_LPP/g' < input_template > input_LPP.temp 
      sed 's/dtTBA/2.0e-7/g' < input_LPP.temp > input_LPP.temp1
      sed 's/sigmaTBA/0.0001/g' < input_LPP.temp1 > input_LPP
      rm input_LPP.temp*
	   ln -s input_LPP input
      sed 's/ForT/T/g' < inputlpp_template > inputlpp_LPP 
	   ln -s inputlpp_LPP inputlpp
      mpirun -np $npstart paris > tmpout_LPP 2>&1
      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout_LPP`
      rm  input inputlpp 
      sed 's/out_dir/out_DNS/g' < input_template > input_DNS.temp 
      sed 's/dtTBA/5.0e-8/g' < input_DNS.temp > input_DNS.temp1
      sed 's/sigmaTBA/0.0001/g' < input_DNS.temp1 > input_DNS
      rm input_DNS.temp*
	   ln -s input_DNS input
      sed 's/ForT/F/g' < inputlpp_template > inputlpp_DNS 
	   ln -s inputlpp_DNS inputlpp
      mpirun -np $npstart paris > tmpout_DNS 2>&1
      echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout_DNS`
fi

sh ./plot.sh
sh ./plot-vel.sh

echo -e "$GREEN" "PASS & Check results in bubble_vel.pdf."  "$NORMAL"


