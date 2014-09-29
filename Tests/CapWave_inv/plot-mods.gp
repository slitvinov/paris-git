set term jpeg size 1200, 1200 font "Arial, 16"
set xrange [0:1]
set yrange [0:1]
set output "mods_par.jpg"
set xlabel " X "
set ylabel " Z "
set grid lw 2
set xtics 0.05
set ytics 0.05
set key right

splot './out/x_mod-00000-0000'.t.'.txt'  ps 2 lc 1 title 'X0', './out/x_mod-00001-0000'.t.'.txt'  ps 2 lc 2 title 'X1', './out/x_mod-00002-0000'.t.'.txt'  ps 2 lc 3 title 'X2', './out/x_mod-00003-0000'.t.'.txt'  ps 2 lc 4 title 'X3', './out/z_mod-00000-0000'.t.'.txt'  pt 10 ps 2 lc 1 title 'Z0', './out/z_mod-00001-0000'.t.'.txt'  pt 6 ps 2 lc 2 title 'Z1', './out/z_mod-00002-0000'.t.'.txt'  ps 2 lc 3 title 'Z2', './out/z_mod-00003-0000'.t.'.txt'  ps 2 lc 4 title 'Z3', './out_serial/x_mod-00000-0000'.t.'.txt'  ps 1 title 'serial x', './out_serial/z_mod-00000-0000'.t.'.txt'  ps 1 title 'serial z'

splot './out/x_mod-00000-0000'.t.'.txt'  ps 2 lc 1 title 'X0', './out/x_mod-00001-0000'.t.'.txt'  ps 2 lc 2 title 'X1', './out/x_mod-00002-0000'.t.'.txt'  ps 2 lc 3 title 'X2', './out/x_mod-00003-0000'.t.'.txt'  ps 2 lc 4 title 'X3', './out_serial/x_mod-00000-0000'.t.'.txt'  ps 1 title 'serial x'

splot './out/z_mod-00000-0000'.t.'.txt'  pt 10 ps 2 lc 1 title 'Z0', './out/z_mod-00001-0000'.t.'.txt'  pt 6 ps 2 lc 2 title 'Z1', './out/z_mod-00002-0000'.t.'.txt'  ps 2 lc 3 title 'Z2', './out/z_mod-00003-0000'.t.'.txt'  ps 2 lc 4 title 'Z3', './out_serial/z_mod-00000-0000'.t.'.txt'  ps 1 title 'serial z'
