#!/bin/bash
#set -x

/bin/rm -fr out input stats
let npstart=4
if [ $# -gt 0 ]; then
    if [ $1 == mono ]; then
	echo "mono"
	let npstart=1
	ln -s testinput.mono input
	precision=0.01
    else
	echo "unknown option" $1
	exit 1
    fi
else
    ln -s testinput input
    precision=0.002
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

awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout

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
	pariscompare output1 Poiseuille_theory $precision
    cd ..
else
    echo "FAIL: directory out not created"
fi

if [ $CGFONTGETGLYPH_PARIS_PROBLEM == 0 ]; then
gnuplot <<EOF &
set size square
set grid
set nolabel
set xrange[-0.1:1.1]
set yrange[-0.05:*]
plot "out/output1" w lp notitle, 0.5*x*(1-x) w l notitle, 0 w l notitle
pause 2
exit
EOF
else
gnuplot <<EOF &
set size square
set grid
set nolabel
set xrange[-0.1:1.1]
set yrange[-0.05:*]
plot "out/output1" w lp notitle, 0.5*x*(1-x) w l notitle, 0 w l notitle
exit
EOF
fi

/bin/rm -f input

