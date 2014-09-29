set term jpeg size 1200, 1200 font "Arial, 16"
set xrange [0:1]
set yrange [0:1]
set output "curve_par.jpg"
set xlabel " X "
set ylabel " Z "
set grid lw 2
set xtics 0.05
set ytics 0.05
set key right

splot './out/kappa-00000-0000'.t.'.txt'  ps 2 lc 1, './out/kappa-00001-0000'.t.'.txt'  ps 2 lc 2, './out/kappa-00002-0000'.t.'.txt'  ps 2 lc 3, './out/kappa-00003-0000'.t.'.txt'  ps 2 lc 4, './out_serial/kappa-00000-0000'.t.'.txt'  ps 2 title 'serial kappa'
