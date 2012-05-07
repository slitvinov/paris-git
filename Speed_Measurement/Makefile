#------------------------------------------------

test1:
	make clean
	ln -s input1 input
	mpirun -np 1 ../run.exe

test2:
	make clean
	ln -s input2 input
	mpirun -np 4 ../run.exe


test3:
	make clean
	ln -s input3 input
	mpirun -np 4 ../run.exe

test4:
	make clean
	ln -s input4 input
	mpirun -np 8 ../run.exe

clean:
	rm -fR out*  input *~ */*~  tmp.* track stats

