set size square
set xtics 0.03125,0.0625,0.96875 format ""
set ytics 1
set xlabel "x"
set ylabel "height in grid units"
set yrange [-8:8]
set grid lt 9
plot "out/reference.txt" w l t "expected values", 0 w l notitle, "out/output1" w l t "calculated values"
set term pdf
set out "heights_test.pdf"
replot

