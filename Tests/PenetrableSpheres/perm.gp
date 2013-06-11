
set grid
set xrange [0:0.98]
set yrange [1e-5:1e2]
set xlabel "solid fraction"
set ylabel "permeability"
#set key left
set log y
plot "permgp.tmp" u (1-$2):3 t  "ParisSimulator" , (1-x)/(54*log(1-x)**2) w lines t "Kozeny-Carman"
