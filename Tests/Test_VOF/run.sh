#!/bin/bash
set -x


echo `pwd`
rm -fR out
mpirun -np 8 paris

