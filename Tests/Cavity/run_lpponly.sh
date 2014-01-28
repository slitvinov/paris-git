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
fi

