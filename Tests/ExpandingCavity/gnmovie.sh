#!/bin/bash
#set -x

LIMITE=20
START=0
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
done

gnuplot << EOF
n=0
call "gnmovie.gp" $LIMITE
call "merged.gp" $LIMITE
exit
EOF

exit 0
