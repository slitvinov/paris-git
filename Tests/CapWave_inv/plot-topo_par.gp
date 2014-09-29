set term jpeg size 1200, 1200 font "Arial, 16"
set xrange [0:1]
set yrange [0:1]
set output "Topo_par.jpg"
set xlabel " X "
set ylabel " Z "
set grid lw 2
set xtics 0.05
set ytics 0.05
unset key

plot './out/Top_1-00000-0000'.t.'.txt' using 1:3 ps 2 lc 1 title 'L1', './out/Top_1-00001-0000'.t.'.txt' using 1:3 ps 2 lc 2, './out/Top_1-00002-0000'.t.'.txt' using 1:3 ps 2 lc 3, './out/Top_1-00003-0000'.t.'.txt' using 1:3 ps 2 lc 4, './out/Top_2-00000-0000'.t.'.txt' using 1:3 ps 2 lc 1, './out/Top_2-00001-0000'.t.'.txt' using 1:3 ps 2 lc 2, './out/Top_2-00002-0000'.t.'.txt' using 1:3 pt 8 ps 2 lc 3, './out/Top_2-00003-0000'.t.'.txt' using 1:3 pt 6 ps 2 lc 4, './out/Top_0-00000-0000'.t.'.txt' using 1:3 pt 6 ps 2 lw 2 lc 1, './out/Top_0-00001-0000'.t.'.txt' using 1:3 ps 2 lc 2, './out/Top_0-00002-0000'.t.'.txt' using 1:3 ps 2 lc 3, './out/Top_0-00003-0000'.t.'.txt' using 1:3 ps 2 lc 4, './out/P1-00000-0000'.t.'.txt' using 1:3 ps 2 pt 7 lc 1, './out/P1-00001-0000'.t.'.txt' using 1:3 ps 2 pt 7 lc 2, './out/P1-00002-0000'.t.'.txt' using 1:3 ps 2 pt 7 lc 3, './out/P1-00003-0000'.t.'.txt' using 1:3 ps 2 pt 7 lc 4
