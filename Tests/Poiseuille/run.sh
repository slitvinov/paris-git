#!/bin/bash
set -x

/bin/rm -fr out input
let npstart=4
ln -s testinput input
if [ `awk 'BEGIN {FS = "=";}  /npx/ {print $2}' < input` == '1' ]; then
  echo "mono"
  let npstart=1
fi

if [ `awk 'BEGIN {FS = "=";}  /DoFront/ {print $2}' < input` == 'T' ]; then
    let np=$npstart+1
else
    let np=$npstart
fi

if [ $np -gt 1 ]; then
    mpirun -np $np paris > tmpout
else
    paris > tmpout
fi

if [ -d out ]; then
    cd out
    if [ $npstart == 4 ]; then
	head -n 3 output_location00000 > output1
	cat output_location00001 >> output1
    elif [ $npstart == 1 ]; then
	cp output_location00000 output1
    else
	echo "case not handled: npstart = ",$npstart
    fi
	compare output1 Poiseuille_theory 0.002
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

