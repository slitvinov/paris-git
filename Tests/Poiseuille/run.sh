#!/bin/bash
#set -x

/bin/rm -fr out input
ln -s testinput input
mpirun -np 4 paris
if [ -d out ]; then
    cd out
    head -n 3 output_location00000 > output1
    cat output_location00001 >> output1
    compare output1 Poiseuille_theory 0.002
    /bin/rm -f output1
    cd ..
else
    echo "FAIL: directory out not created"
fi
sh ./plot.sh
/bin/rm -f input

