set log x
set log y
set xlabel "Grid points per Radius"
set ylabel "Curvature error in 3D"
plot "gerris-3D.tmp" u 1:2 t "Gerris",  "paris-3D-2.tmp" u 1:2 t "ParisSim depth=2 with mixed heights",  "paris-3D-3.tmp" u 1:2 t "ParisSim depth=3 with mixed heights", 2/(x*x) # ,   "paris-3D-2-nomixed.txt" u 1:2 t "ParisSim depth=2 with mixed heights",  "paris-3D-3-nomixed.txt" u 1:2 t "ParisSim depth=3 with mixed heights"
set term pdf
set out "curvature-3D.pdf"
plot "gerris-3D.tmp" u 1:2 t "Gerris",  "paris-3D-2.tmp" u 1:2 t "ParisSim depth=2 with mixed heights",  "paris-3D-3.tmp" u 1:2 t "ParisSim depth=3 with mixed heights", 2/(x*x) 





