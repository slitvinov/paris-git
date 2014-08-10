set xlabel "time"
set ylabel "w(0,0,R)"
set yrange [*:*]
plot "reference.txt" t "reference" w l, "test-vel-NH-NI-notMOF-lotol.txt" t "low tolerance VOF" w l, "out/droplet-test-vel.txt" t "test simulation nonMOF" w l, 0 notitle
#plot "out/droplet-test-vel.txt" t "with Hypre" w l, "test-vel-NH-NI-notMOF-lotol.txt" t "low tolerance" w l, 0 notitle
set term pdf
set out "tmp.pdf"
replot


