set log x
set log y
set xlabel "Grid points per Radius"
set ylabel "Curvature error in 3D"
plot "gerris-3D.tmp" u 1:2 t "Gerris",  "paris-3D.tmp" u 1:2 t "ParisSimulator", 2/(x*x) 
set term pdf
set out "curvature-3D.pdf"
plot "gerris-3D.tmp" u 1:2 t "Gerris",  "paris-3D.tmp" u 1:2 t "ParisSimulator", 2/(x*x) 





