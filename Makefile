#--------- paris-ParisSimulator main Makefile --------------------------

# babbage
OMPI_FC=gfortran44

FC = mpif90

# remove funny cflags from my environment

FFLAGS = -O3 # -g -gstabs
CFLAGS = -g -gstabs

# select option for hypre
# default hypre installation without root privileges:
# HYPRE_DIR = $(HOME)/hypre-2.8.0b/src
# Macport installation in /opt 
HYPRE_DIR = /opt/hypre
# babbage
# HYPRE_DIR = /share/apps/hypre

HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 


#------------------------No changes needed beyond this line----------------------------------------------


OBJ = paris.o solids.o modules.o vofmodules.o front.o
OBJC = compare.o
SRC = $(wildcard  *.f90) 

install: $(OBJ)
	@echo compiler is FC = $(FC), mpi override is OMPI_FC = $(OMPI_FC)
	$(FC) -o paris $(FOPTS) $(OBJ) $(FOBJ) $(HYPRE_LIBS) 
	@if [ ! -d $(BINDIR) ] ; then echo "directory bin does not exist creating it" ; mkdir $(BINDIR) ; fi 
	mv paris $(BINDIR)/paris	

all: tags install compare

clean:
	@rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.*
	@cd Tests; sh ./clean.sh; cd ..
	@cd Documentation; make clean; cd ..

distclean: clean
	@rm -fR  session* *.xml TAGS tags input

test:  install compare
	@cd Tests; sh ./runtests.sh

# single processor test
minitest: install
	cd Tests/Mini; sh run.sh

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 

paris.o:  paris.f90 solids.o modules.o vofmodules.o front.o
	$(FC) -c  $(FFLAGS) $<

vofmodules.o: vofmodules.f90 modules.o
	$(FC) -c  $(FFLAGS) $<

solids.o:  solids.f90 modules.o
	$(FC) -c   $(FFLAGS) $<

front.o:  front.f90 modules.o
	$(FC) -c   $(FFLAGS) $<

%.o : %.f90
	$(FC) -c  $(FFLAGS) $<

compare: $(OBJC)
	$(CC) -o compare  $(OBJC) -lm
	mv compare ~/bin

.c.o:   $< 
	cc -c $(CFLAGS)   $< 
