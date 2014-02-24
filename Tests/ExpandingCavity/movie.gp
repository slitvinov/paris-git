set term gif
set output 'bubble.gif'

set xrange [0:1];
set yrange [0:1];

set key below;
set size square;

set xlabel " x "
set ylabel " y "

unset surface
set contour
set cntrparam linear
set view map
set cntrparam levels discrete 0.75

splot "toplot.tmp" w l notitle
