unset log y
set log x

set yrange [0:1]
set ylabel "occurence of method"
set xlabel "Grid points per Radius"
plot "method_count-2.tmp" u 1:2 t "9 heights depth=2" w lp,  "method_count-3.tmp" u 1:2 t "9 heights depth=3" w lp,
set term pdf
set out "curvature-method.pdf"
# plot "method_count.tmp" u 1:2 t "9 heights" w lp, "method_count.tmp" u 1:3 t "mixed heights" w lp,  "method_count.tmp" u 1:4 t "centroids" w lp
plot "method_count-2.tmp" u 1:2 t "9 heights depth=2" w lp,  "method_count-3.tmp" u 1:2 t "9 heights depth=3" w lp,
