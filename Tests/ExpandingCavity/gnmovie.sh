#!/bin/bash
#set -x

LIMITE=20
START=0
DELTA=1000

/bin/rm -f *.png

for ((a=START; a <= LIMITE ; a++)) # Doubles parenthèses, et "LIMITE" sans "$".
do
#let "p=$p+1"
#echo $a
let frame=$DELTA*$a
if [ $a == 0 ]; then
    frame=0000
fi
if [ $a -lt 10 ]; then
    frame='0'$frame
fi

#echo "frame" $frame
cp "out/CVoF-00000-$frame.txt" toplot$a.tmp
done

gnuplot <<EOF
n=0
call "gnmovie.gp" $a
exit
EOF

exit 0
