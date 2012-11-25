
mono64-1:
	make clean
	ln -s input-mono64-1 input
#	$(OPENMPIV)/bin/mpirun -np 1 $(HOME)/bin/paris > out-mono64-1 2>err-mono64-1
	mpirun -np 1 $(HOME)/bin/paris > out-mono64-1 2>err-mono64-1


mono64-8:
	make clean
	ln -s input-mono64-8 input
#	$(OPENMPIV)/bin/mpirun -np 8 $(HOME)/bin/paris > out-mono64-8 2>err-mono64-8
	mpirun -np 8 $(HOME)/bin/paris > out-mono64-8 2>err-mono64-8


mono64-64:
	make clean
	ln -s input-mono64-64 input
#	$(OPENMPIV)/bin/mpirun -np 64 $(HOME)/bin/paris > out-mono64-64 2>err-mono64-64
	mpirun -np 64 $(HOME)/bin/paris > out-mono64-64 2>err-mono64-64


mono64-512:
	make clean
	ln -s input-mono64-512 input
#	$(OPENMPIV)/bin/mpirun -np 512 $(HOME)/bin/paris > out-mono64-512 2>err-mono64-512
	mpirun -np 512 $(HOME)/bin/paris > out-mono64-512 2>err-mono64-512


mono64-4096:
	make clean
	ln -s input-mono64-4096 input
#	$(OPENMPIV)/bin/mpirun -np 4096 $(HOME)/bin/paris > out-mono64-4096 2>err-mono64-4096
	mpirun -np 4096 $(HOME)/bin/paris > out-mono64-4096 2>err-mono64-4096


mono64-32768:
	make clean
	ln -s input-mono64-32768 input
#	$(OPENMPIV)/bin/mpirun -np 32768 $(HOME)/bin/paris > out-mono64-32768 2>err-mono64-32768
	mpirun -np 32768 $(HOME)/bin/paris > out-mono64-32768 2>err-mono64-32768

#------------------

clean:
	rm -fr out* err* stats track input *~ #*

