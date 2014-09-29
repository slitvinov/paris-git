set term jpeg size 1200, 1200 font "Arial, 16"
set xrange [0:1]
set yrange [0:1]
set output "phase_par.jpg"
set xlabel " X "
set ylabel " Z "
set grid lw 2
set xtics 0.05
set ytics 0.05
set key right

splot './out/phase-00000-0000'.t.'.txt'  ps 2 lc 1 title 'X0', './out/phase-00001-0000'.t.'.txt'  ps 2 lc 2 title 'X1', './out/phase-00002-0000'.t.'.txt'  ps 2 lc 3 title 'X2', './out/phase-00003-0000'.t.'.txt'  ps 2 lc 4 title 'X3', './out_serial/phase-00000-0000'.t.'.txt'  ps 2 title 'serial x', './out_serial/z_mod-00000-0000'.t.'.txt'  ps 2 title 'serial z'
