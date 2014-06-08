set size square
set nolabel
set xrange[*:*]
set yrange[*:*]


set xlabel " i "
set ylabel " j "

unset surface
set contour
set cntrparam linear
set view map
set cntrparam levels auto 10

set title "pressure"
# set key off

# plot "2plot.txt" w vec notitle
splot "presplot.txt" w l 

exit


