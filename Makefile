#--------- paris-ParisSimulator main Makefile --------------------------

# babbage
OMPI_FC=gfortran44

FC = mpif90
CC = mpicc

# remove funny cflags from my environment
CFLAGS = -O
FFLAGS = -O3

# select option for hypre
# default hypre installation without root privileges:
# HYPRE_DIR = $(HOME)/hypre-2.8.0b/src
# Macport installation in /opt 
HYPRE_DIR = /opt/hypre
# babbage
# HYPRE_DIR = /share/apps/hypre

HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 


#------------------------No changes needed beyond this line----------------------------------------------


OBJ = paris.o utilc.o solids.o modules.o vofmodules.o
FOBJ = # flux3d.o alpha3d.o
SRC = $(wildcard *.h *.c *.f90) 
INC = $(wildcard *.h) 

install: $(OBJ)
	@echo compiler is FC = $(FC), mpi override is OMPI_FC = $(OMPI_FC)
	$(FC) -g -o paris $(FOPTS) $(OBJ) $(FOBJ) $(HYPRE_LIBS) 
	mv paris ~/bin/paris

all: tags install

clean:
	@rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.*

distclean:
	@rm -fR *.o *.mod paris *.gz stats *~ track out* errftc tmp* *.tmp fort.* *.visit TAGS tags core.* input
	@cd Speed_Measurement; make clean; cd ..
	@cd Poiseuille_Test; make clean; cd ..
	@cd Test_VOF; make clean; cd ..
	@cd Documentation; make clean; cd ..

test:	install
	@rm -fR out input
	@ln -s miniinput input
	mpirun -np 8 paris

# single processor test
minitest:	install
	@rm -fR out input
	@ln -s minimono input
	paris

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 

paris.o:  paris.f90 solids.o modules.o vofmodules.o $(FOBJ)
	$(FC) -c  $(FFLAGS) $<

vofmodules.o: vofmodules.f90 modules.o
	$(FC) -c  $(FFLAGS) $<

solids.o:  solids.f90 modules.o
	$(FC) -c   $(FFLAGS) $<

%.o : %.f90
	$(FC) -c  $(FFLAGS) $<


.c.o:   $< $(INC)
	$(CC) -c $(CFLAGS)   $< 


