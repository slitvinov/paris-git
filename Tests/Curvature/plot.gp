set log x
set log y
set xlabel "Grid points per Radius"
set ylabel "Curvature error in 3D"
plot "paris-3D-3.tmp" u 1:2 t "ParisSimulator", 2/(x*x) axes x1y1
#set term pdf
#set multiplot
#set out "curvature-3D.pdf"
#set yrange [*:*] 
#set log y
#set xlabel "Grid points per Radius"
#set ylabel "Curvature error in 3D"
#plot "paris-3D.tmp" u 1:2 t "ParisSimulator", 2/(x*x)
#unset log y
#set yrange [0:1]
#plot "method_count.tmp" u 1:2 t "9 heights", "method_count.tmp" u 1:3 t "mixed heights", "method_count.tmp" u 1:4 t "centroids"





