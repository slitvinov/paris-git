set term jpeg size 1200, 1200 font "Arial, 16"
set xrange [0:1]
set yrange [0:1]
set output "norms_par.jpg"
set xlabel " X "
set ylabel " Z "
set grid lw 2
set xtics 0.05
set ytics 0.05
set key right

splot './out/n1-00000-0000'.t.'.txt'  ps 2 lc 1 title 'n1-0', './out/n1-00001-0000'.t.'.txt'  ps 2 lc 2 title 'n1-1', './out/n1-00002-0000'.t.'.txt'  ps 2 lc 3 title 'n1-2', './out/n1-00003-0000'.t.'.txt'  ps 2 lc 4 title 'n1-3', './out/n3-00000-0000'.t.'.txt'  pt 12 ps 2 lc 1 title 'n3-0', './out/n3-00001-0000'.t.'.txt'  ps 2 lc 2 title 'n3-1', './out/n3-00002-0000'.t.'.txt'  pt 8 ps 2 lc 3 title 'n3-2', './out/n3-00003-0000'.t.'.txt'  pt 6 ps 2 lc 4 title 'n3-3', './out_serial/n1-00000-0000'.t.'.txt'  ps 2 title 'serial n1', './out_serial/n3-00000-0000'.t.'.txt'  ps 2 title 'serial n3'
