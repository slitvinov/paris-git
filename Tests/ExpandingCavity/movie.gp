set term png
set output 'bubble.png'

set xrange [0:1]
set yrange [0:1]

#set key 2,3.5,2
set key below
set size square
#set  origin 0,0.66

set xlabel " x "
set ylabel " y "

unset surface
set contour
set cntrparam linear
set view map
set cntrparam levels discrete 0.5

splot "toplot.tmp" w l notitle
