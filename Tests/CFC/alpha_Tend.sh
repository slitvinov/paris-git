#! /bin/bash
set -x

endtime=`awk '/EndTime/ { print $3 }' < input`
echo alpha Tend
awk -v  endtime=$endtime '{print  endtime/$2}' < out/flowrate.txt
cat > tmp.txt <<EOF  
  9.322244E-05  9.322372E-05   9.322498E-05  16 1 0.1e-3
  8.997408E-05  8.997437E-05   8.997466E-05 32 1 0.1e-3 
  8.528072E-05  8.527347E-05  8.526733E-05  64  5 0.04e-3
EOF
   
echo alpha Tend from stats 
awk  -v endtime=$endtime '{ print $4 " : " ((1. -  ($3-$2)/($2 - $1))/($5*$6)) }' < tmp.txt
#awk -v dt=$dt -v endtime=$endtime '{ print ($2 - $1)/($3-$2) }' < tmp.txt
  # exp ( -alpha delta t)  = 1 - alpha delta t
