set log x
set log y
set xlabel "Grid points per Radius"
set ylabel "Curvature error in 3D"
plot "gerris-3D.tmp" u 1:2 t "Gerris L2 one case",  "paris-3.tmp" u 1:2 t "ParisSim average L2", 1/(x*x) t "x^-2",  "paris-3.tmp" u 1:3 t "ParisSim max Linf", 1/x t "x^-1"
set term pdf
set out "curvature-3D-inf.pdf"
plot "gerris-3D.tmp" u 1:2 t "Gerris L2 one case",  "paris-3.tmp" u 1:2 t "ParisSim average L2", 1/(x*x) t "x^-2",  "paris-3.tmp" u 1:3 t "ParisSim max Linf", 1/x t "x^-1"





