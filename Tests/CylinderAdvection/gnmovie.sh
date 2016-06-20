#! /bin/bash
#set -x

LIMITE=38
START=1
DELTA=100

/bin/rm -f *.png

for ((a=START; a <= LIMITE ; a++)) # Doubles parenthèses, et "LIMITE" sans "$".
do

let frame=$DELTA*$a
if [ $a == 0 ]; then
    frame=000
fi
if [ $a -lt 10 ]; then
    frame='00'$frame
else
    frame='0'$frame
fi

cp "out/CVoF-00000-$frame.txt" toplot$a.tmp
cp toplot$a.tmp toplot.tmp
gnuplot < gnframe.gp
cp tmp.png frame$frame.png
echo -n "."
done
echo " "
convert frame*.png  CAmovie.gif

exit 0
