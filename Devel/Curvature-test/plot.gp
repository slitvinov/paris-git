# may not work because of the shell variables.
set log x
set log y
set xlabel "Grid points per Diameter"
set ylabel "Curvature error in 3D"
plot "gerris-$dim.txt" u (2*$1):2 t "Gerris L2 one case",  "paris-$dim-$ndepth.tmp" u (2*$1):2 t "ParisSim average L2", 4/(x*x) t "x^-2",  "paris-$dim-$ndepth.tmp" u (2*$1):3 t "ParisSim max Linf", 2/x t "x^-1"
set term pdf
set out "curvature-3D.pdf"
replot





