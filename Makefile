#------------------------------------------------
#-------------- Uses openmpi library --------------
#------------------------------------------------

FFLAGS   = -O3 
FC   = mpif90
CFLAGS = -O -I/opt/local/include/openmpi/

# default hypre installation without root privileges:
# HYPRE_DIR = $(HOME)/hypre-2.8.0b/src
# installation in /opt with root
HYPRE_DIR = /opt/hypre
HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 
# babbage
# HYPRE_DIR = /share/apps/hypre

#-----------------------------------------------------------------------------------------------


OBJ     = ftc3d2011.o utilc.o solids.o modules.o
SRC = $(wildcard *.h *.c *.f90) 
INC = $(wildcard *.h) 

install: $(OBJ)
	$(FC) -o run.exe $(FOPTS) $(OBJ) $(HYPRE_LIBS) 
#	put it in a sound location @SZ
	mv run.exe ~/bin/ftc3d2011

all: tags install

clean:
	rm -fR *.o *.mod run.exe *.gz stats *~ track out* stats errftc outftc tmp *.visit TAGS tags
	cd Speed_Measurement; make clean; cd ..
	cd Poiseuille_Test; make clean; cd ..

test:	install
	rm -fR out input
	ln -s miniinput input
	mpirun -np 8 ftc3d2011

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags ftc3d2011.f90 

ftc3d2011.o:  ftc3d2011.f90 solids.o modules.o
	$(FC) -c $(FFLAGS) ftc3d2011.f90 

solids.o:  solids.f90 modules.o
	$(FC) -c $(FFLAGS) $<


%.o : %.f90
	$(FC) -c $(FFLAGS) $<

.c.o:   $< $(INC)
	$(CC) -c $(CFLAGS)   $< 


