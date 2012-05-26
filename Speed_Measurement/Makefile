##------------------------------------------------
test64:
	make clean
	ln -s input64 input
	$(OPENMPIV)/bin/mpirun -np 64 $(HOME)/bin/ftc3d2011 > outftc 2>errftc


test1:
	make clean
	ln -s input1 input
	mpirun -np 1 ftc3d2011

test2:
	make clean
	ln -s input2 input
	mpirun -np 8 ftc3d2011


test-np4:
	make clean
	ln -s input-np4 input
	mpirun -np 4 ftc3d2011


test3:
	make clean
	ln -s input3 input
	mpirun -np 32 ftc3d2011

#test4:
test32:
	make clean
	ln -s input4 input
	$(OPENMPIV)/bin/mpirun -np 32 $(HOME)/bin/ftc3d2011 > outftc 2>errftc

test5:
	make clean
	ln -s input5 input
	mpirun -np 16 ftc3d2011


test6:
	make clean
	ln -s input6 input
	mpirun -np 4 ftc3d2011

clean:
	rm -fR out*  input *~ */*~  tmp.* track stats

