#! /bin/bash
#set -x

if [ $# -lt 1 ]; then
    echo "$0: missing arguments"
    echo usage $0 FILE
    exit
fi

# For CFC
#nmeshes=8
#size=1e-4
#radius=`awk -v size=$size -v nmeshes=$nmeshes 'BEGIN { meshsize=size/nmeshes; print (sqrt(2.)/4)*meshsize}'`


#rescale box to 1
scale=1e4
# radius=6.25e-6
radius=6.25e-6
shift=0.5d0

awk 'END { print "NumSpheres = " $1"\n" }' < $1 > inputcenters.tmp

awk -v shift=$shift -v radius=$radius -v scale=$scale '{ print "sxyzrad(1,"$1") = "$2*scale+shift "\nsxyzrad(2,"$1") = "$3*scale+shift "\nsxyzrad(3,"$1") = "$4*scale+shift "\nsxyzrad(4,"$1") = "radius*scale }' < $1 >> inputcenters.tmp