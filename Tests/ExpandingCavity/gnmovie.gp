#set term png
#set output 'bubble.png'

#set xrange [0:1]
#set yrange [0:1]

set size square

set xlabel " i "
set ylabel " j "

unset surface
set contour
set cntrparam linear
set view map
set cntrparam levels discrete 0.5
set title "ouptut nr  ".n

splot "toplot".n.".tmp" w l notitle lt n

#print  "toplot".n.".tmp  plotted"
pause 1
n=n+1

if (n < 20 ) reread
exit


