#set term png
#set output 'bubble.png'

#set xrange [0:1]
#set yrange [0:1]

set size square

set xlabel " i "
set ylabel " j "
set term pngcairo
set out "tmp.png"
unset surface
set contour
set cntrparam linear
set view map
set cntrparam levels discrete 0.5
#set title "ouptut nr  ".frame

splot "toplot.tmp" w l notitle lt 1
exit


