#--------- paris-ParisSimulator main Makefile --------------------------

# babbage has several fortran compilers
OMPI_FC=gfortran44

FC = mpif90
CC = mpicc

# remove funny cflags from gerris in my environment
CFLAGS = 

# select option for hypre
# default hypre installation without root privileges:
# HYPRE_DIR = $(HOME)/hypre-2.8.0b/src
# Macport installation in /opt 
HYPRE_DIR = /opt/hypre
# babbage
# HYPRE_DIR = /share/apps/hypre

#------------------------ No changes needed beyond this line --------------------------------

HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 

OBJ = paris.o utilc.o solids.o modules.o
SRC = $(wildcard *.h *.c *.f90) 
INC = $(wildcard *.h) 

install: $(OBJ)
	@echo compiler is FC = $(FC), mpi override is OMPI_FC = $(OMPI_FC)
	$(FC) -o paris $(FOPTS) $(OBJ) $(HYPRE_LIBS) 
	mv paris ~/bin/paris

all: tags install

clean:
	@rm -fR *.o *.mod paris *.gz stats *~ track out* stats errftc tmp* *.tmp fort.* *.visit TAGS tags core.*
	@cd Speed_Measurement; make clean; cd ..
	@cd Poiseuille_Test; make clean; cd ..

test:	install
	@rm -fR out input
	@ln -s miniinput input
	mpirun -np 8 paris

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 

paris.o:  paris.f90 solids.o modules.o
	$(FC) -c  $<

solids.o:  solids.f90 modules.o
	$(FC) -c  $<

%.o : %.f90
	$(FC) -c $<


.c.o:   $< $(INC)
	$(CC) -c $(CFLAGS)   $< 


