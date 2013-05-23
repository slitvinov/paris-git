#! /bin/bash
set -x

endtime=`awk '/EndTime/ { print $3 }' < input`
echo alpha Tend
awk -v  endtime=$endtime '{print  endtime/$2}' < out/flowrate.txt
cat > tmp.txt <<EOF  
 8.528072E-05  8.527347E-05  8.526733E-05  
EOF
   
echo alpha Tend from stats 
dt=0.2e-3
awk -v dt=$dt -v endtime=$endtime '{ print ((1. -  ($3-$2)/($2 - $1))/dt)*endtime }' < tmp.txt
#awk -v dt=$dt -v endtime=$endtime '{ print ($2 - $1)/($3-$2) }' < tmp.txt
  # exp ( -alpha delta t)  = 1 - alpha delta t
