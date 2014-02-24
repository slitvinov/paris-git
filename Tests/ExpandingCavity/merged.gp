set term gif
set output 'bubble.gif'

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
lim=$0
set title "ouptut from 0 to ".lim
set key off

splot for [n=1:lim] "toplot".n.".tmp" w l notitle 

exit


