set xlabel "time"
set ylabel "w(0,0,5R/4)"
set yrange [*:*]
plot "test-vel-NH-NI-notMOF-lotol.txt" t "Nonmixed" w l, 0 notitle, "out/droplet-test-vel.txt" t "mixed" w l
#plot "reference.txt" t "reference" w l, "test-vel-NH-NI-notMOF-lotol.txt" t "low tolerance VOF" w l, "out/droplet-test-vel.txt" t "test simulation nonMOF" w l, 0 notitle
#plot "out/droplet-test-vel.txt" t "with Hypre" w l, "test-vel-NH-NI-notMOF-lotol.txt" t "low tolerance" w l, 0 notitle
set term pdf
set out "tmp.pdf"
replot


# N=32 R=0.3 K=0.75   3/4 / 3/5 = 5/4