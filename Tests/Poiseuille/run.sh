#!/bin/bash
set -x

/bin/rm -fr out input
let npstart=4
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

gnuplot <<EOF
set size square
set grid
set nolabel
set xrange[-0.1:1.1]
set yrange[-0.05:0.15]
plot "out/output1" w lp notitle, 0.5*x*(1-x) w l notitle, 0 w l notitle
set term pdf
set out 'Poiseuille_plot.pdf'
plot "out/output1" w lp notitle, 0.5*x*(1-x) w l notitle, 0 w l notitle
exit
EOF

/bin/rm -f input

