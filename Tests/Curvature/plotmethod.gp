set log y
set log x

set title "Method distribution for depth of stack NDEPTH = ndepth"
set xrange [1:110]
set yrange [0.000001:2]
set ylabel "occurence of centroid method"
#set xlabel "Grid points per Radius"
#plot "method_count-3D-3.tmp" u ($1):3 t "centroids" , 1

set xlabel "Grid points per Diameter"
#plot "method_count-3D-3.tmp" u (2*$1):3 t "centroids" , 1 t "100 %"

# plot "method_count.tmp" u 1:2 t "9 heights" w lp
#plot "method_count-3D-3.tmp" u 1:2 t "9 heights", "method_count-3D-3.tmp" u 1:3 t "mixed heights",  "method_count-3D-3.tmp" u 1:4 t "centroids"
plot "method_count-3D-ndepth.tmp" u 1:2 t "9 heights", "method_count-3D-ndepth.tmp" u 1:3 t "mixed heights",  "method_count-3D-ndepth.tmp" u 1:4 t "centroids"
# plot "method_count-2.tmp" u 1:2 t "9 heights depth=2" w lp,  "method_count-3.tmp" u 1:2 t "9 heights depth=3" w lp,
set term pdf
set out "curvature-method.pdf"
replot
