test:
	make clean
	ln -s testinput input
	mpirun -np 4 ftc3D2011
	sh ./plot.sh

clean:
	rm -fR out* input *~  track stats tmp.* *.visit core.* \#*

