#!/bin/bash
#set -x
rm -fr cav_vol*
./run-3dtest.sh 1 1 1 30 #this arrangement must be here in order to have a serial reference
./run-3dtest.sh 3 1 1 30
./run-3dtest.sh 4 1 1 32
./run-3dtest.sh 2 1 4 32
#./run-3dtest.sh 3 3 3 30

#$1 = NPX number of procs in x-direction
#$2 = NPY
#$3 = NPZ
#$4 = NX grid points per coord direction
