#! /bin/bash
#set -x

if [ $# -lt 1 ] 
then
 echo "usage   : " $0 "itimestep  [processor rank]"
 echo "example : " $0 "987 0"
 exit
fi
echo arg 1 is $1

if  [ $# -lt 2 ] 
then
  rm -f 2plot.txt
  padding="0000"
  if [ 2 -eq `echo '$1'| wc -c` ] 
      then
      padding="000"
  elif [ 3 -eq `echo '$1'| wc -c` ]; then
      padding="00"
  elif [ 4 -eq `echo '$1'| wc -c` ]; then
      padding="0"
  elif [ 5 -eq `echo '$1'| wc -c` ]; then
      padding=""
  else
      echo argument 1 too long
      exit 1
  fi
      cat  out/UV-00000-$padding$1.txt >> 2plot.txt   
      cat  out/UV-00001-$padding$1.txt >> 2plot.txt   
      cat  out/UV-00002-$padding$1.txt >> 2plot.txt   
      cat  out/UV-00003-$padding$1.txt >> 2plot.txt   
else
    cp out/UV-0000$2-$padding$1.txt 2plot.txt
fi
gnuplot <<EOF
set size square
set nolabel
#set xrange[-0.5:0.5]
#set yrange[-0.5:0.5]
plot "2plot.txt" w vec notitle
EOF