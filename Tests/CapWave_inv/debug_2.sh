#!/bin/bash

gnuplot <<EOF
	t=1
	call "plot-phase.gp"
	call "plot-norm.gp"
	call "plot-curve.gp"
EOF
