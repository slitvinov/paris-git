set size square
set xtics 0.03125,0.0625,0.96875 format ""
set ytics 1
set xlabel "x"
set ylabel "height in grid units"
set yrange [-8:8]
set grid lt 9
plot "out/reference.txt" t "expected values", 0 w l notitle, "out/output1a" w l t "calculated values in proc/subd above", "out/output1b" w l t "calculated values in proc/subd below"

