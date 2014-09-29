#!/bin/bash

gnuplot <<EOF
	t=1
	#call "plot-topo_par.gp"
	call "plot-pg.gp"
	call "plot-mods.gp"
EOF

