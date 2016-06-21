#set term png
#set output 'bubble.png'

#set xrange [0:1]
#set yrange [0:1]

set size square

set xlabel " x "
set ylabel " y "
set term pngcairo
set out "tmp.png"
unset surface
set contour
set cntrparam linear
set view map
unset clabel
set cntrparam levels discrete 0.5

#set title "ouptut nr  ".frame
set grid

splot "toplot.tmp" w l  notitle 
exit


