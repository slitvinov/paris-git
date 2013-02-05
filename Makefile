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
	@cd Tests; sh ./clean.sh; cd ..
	@cd Documentation; make clean; cd ..

distclean: clean
	@rm -fR  session* *.xml TAGS tags input

test:  install
	@cd Tests; sh ./runtests.sh

# single processor test
minitest: install
	cd Tests/Mini; rm -fR out input; ln -s minimono input; paris

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


