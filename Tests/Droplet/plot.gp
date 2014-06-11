set xlabel "time"
set ylabel "w(0,0,R)"
set yrange [-0.0002:0.0002]

set term pdf
set out "tmp.pdf"
#plot "out/droplet-test-vel.txt" t "test simulation" w l, "reference.txt" t "reference" w l, 0 notitle
plot "out/droplet-test-vel.txt" t "with Hypre" w l, "test-vel-NH-NI-notMOF-lotol.txt" t "low tolerance" w l, 0 notitle
set term aqua
replot


