
set grid
set xrange [0:0.98]
set yrange [1e-5:1e2]
set xlabel "solid fraction"
set ylabel "permeability K/R^2"
#set key left
set log y
#plot "permgp.tmp" u (1-$2):4 t  "ParisSimulator phi1" ,  "permgp.tmp" u (1-$3):4 t  "ParisSimulator phi2" , (1-x)/(54*log(1-x)**2) w lines t "Kozeny-Carman"
plot "permgp.tmp" u (1-$2):3 t  "ParisSimulator phi1" ,  (1-x)/(54*log(1-x)**2) w lines t "Kozeny-Carman"
