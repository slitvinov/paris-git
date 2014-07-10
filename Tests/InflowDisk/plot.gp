set yrange[1.51:1.65]
set grid
set xlabel "time"
set ylabel "pressure gradient"
plot "stats" u 1:($18 - $19) w l 
set term pdf
set out 'tmp.pdf'
replot
exit
