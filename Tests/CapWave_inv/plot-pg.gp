set term jpeg size 1200, 1200 font "Arial, 16"
set xrange [0:1]
set yrange [0:1]
set output "Pg_par.jpg"
set xlabel " X "
set ylabel " Z "
set grid lw 2
set xtics 0.05
set ytics 0.05
set key right

splot './out/P_int1-00000-0000'.t.'.txt'  ps 2 lc 1 title 'P1', './out/P_int1-00001-0000'.t.'.txt'  ps 2 lc 2 title 'P1', './out/P_int1-00002-0000'.t.'.txt'  ps 2 lc 3 title 'P1', './out/P_int1-00003-0000'.t.'.txt'  ps 2 lc 4 title 'P1', './out/P_int3-00000-0000'.t.'.txt' pt 7 ps 2 lc 1 title 'P3', './out/P_int3-00001-0000'.t.'.txt'  ps 2 lc 2 title 'P3', './out/P_int3-00002-0000'.t.'.txt'  pt 8 ps 2 lc 3 title 'P3', './out/P_int3-00003-0000'.t.'.txt'  pt 6 ps 2 lc 4 title 'P3', './out_serial/P_int1-00000-0000'.t.'.txt' ps 1 title 'serial Px', './out_serial/P_int3-00000-0000'.t.'.txt' ps 1 lc 0 title 'serial Pz'

splot './out/P_int1-00000-0000'.t.'.txt'  ps 2 lc 1 title 'P1', './out/P_int1-00001-0000'.t.'.txt'  ps 2 lc 2 title 'P1', './out/P_int1-00002-0000'.t.'.txt'  ps 2 lc 3 title 'P1', './out/P_int1-00003-0000'.t.'.txt'  ps 2 lc 4 title 'P1', './out_serial/P_int1-00000-0000'.t.'.txt' ps 1 lc 1 title 'serial Px'

splot './out/P_int3-00000-0000'.t.'.txt' pt 7 ps 2 lc 1 title 'P3', './out/P_int3-00001-0000'.t.'.txt'  ps 2 lc 2 title 'P3', './out/P_int3-00002-0000'.t.'.txt'  pt 8 ps 2 lc 3 title 'P3', './out/P_int3-00003-0000'.t.'.txt'  pt 6 ps 2 lc 4 title 'P3', './out_serial/P_int3-00000-0000'.t.'.txt' ps 1 lc 0 title 'serial Pz'
