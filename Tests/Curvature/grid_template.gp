#set term x11
set term aqua size 846,646
set xrange [0.:0.8]
set yrange [0.:0.8]
set xrange [0.5:0.64]
set yrange [0.5937:0.5938]
set xrange [0.5937:0.5938]
set xrange [0.56:0.64]
set yrange [0.56:0.63]
set multiplot
plot "out/grid.txt" w vec nohead lt 0 notitle
set parametric
set trange [0:2*pi]
# Parametric functions for a circle
fx(t) = RADIUS*cos(t) + XC1
fy(t) = RADIUS*sin(t) + XC2
plot "out/segments.txt" w vec nohead lt 3 notitle
plot  fx(t),fy(t) lt 1 t "exact", "out/points.txt" lw 2 lt 2 t "centroids","out/xheights.txt" pt 10 lt 3 t "x-heights","out/yheights.txt" pt 12 lt 4 t "y-heights"

# pause 100
