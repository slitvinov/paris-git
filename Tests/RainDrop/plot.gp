set xlabel "time"
set ylabel "dimensionless liquid kinetic energy"
set yrange [*:*]
set xrange [0:0.02]
plot "reference-E-k2.txt" t "reference" w l, "E-k2.tmp" t "E-k2" w p, 0 notitle
set term pngcairo
set out "raindrop-E-k2.png"
replot
