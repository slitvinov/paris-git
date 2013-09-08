#--------- paris-ParisSimulator main Makefile --------------------------

# babbage
# OMPI_FC=gfortran44

FC = mpif90

# remove funny cflags from my environment

FFLAGS =  -g -gstabs # -O3 #
CFLAGS = -O # -g -gstabs
BINDIR = $(HOME)/bin

# select option for hypre
# default hypre installation without root privileges:
# HYPRE_DIR = $(HOME)/hypre-2.8.0b/src
# Macport installation in /opt 
HYPRE_DIR = /opt/hypre
# babbage
# HYPRE_DIR = /share/apps/hypre

HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 


#------------------------No changes needed beyond this line----------------------------------------------


OBJ = paris.o solids.o modules.o vofmodules.o front.o surface_tension.o lag_particle.o # uzawa.o
SRC = $(wildcard  *.f90) 

install: $(OBJ)
#	@echo compiler is FC = $(FC), mpi override is OMPI_FC = $(OMPI_FC)
	$(FC) -g -o paris $(FOPTS) $(OBJ) $(FOBJ) $(HYPRE_LIBS) 
	@if [ ! -d $(BINDIR) ] ; then echo "directory bin does not exist creating it" ; mkdir $(BINDIR) ; fi 
	mv paris $(BINDIR)/paris
	@find .  -name "*.sh" -exec chmod +x  {} \; 

all: tags install compare parisdeconv

clean:
	@rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.* statsbub
	@cd Tests; sh ./clean.sh; cd ..
	@cd Documentation; make clean; cd ..

distclean: clean
	@rm -fR  session* *.xml TAGS tags input

test:  install compare parisdeconv
	@echo "This test takes approximately 1 minute on a 4-core intel i7 MacBookPro"
	@cd Tests; chmod +x ./runtests.sh; ./runtests.sh

# single processor test
minitest: install
	@cd Tests/Mini; chmod +x ./run.sh; ./run.sh

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 

paris.o:  paris.f90 solids.o modules.o vofmodules.o front.o surface_tension.o lag_particle.o
	$(FC) -c  $(FFLAGS) $<

vofmodules.o: vofmodules.f90 modules.o
	$(FC) -c  $(FFLAGS) $<

lag_particle.o: lag_particle.f90 vofmodules.o modules.o
	$(FC) -c  $(FFLAGS) $<

surface_tension.o: surface_tension.f90 vofmodules.o modules.o
	$(FC) -c  $(FFLAGS) $<

solids.o:  solids.f90 modules.o
	$(FC) -c   $(FFLAGS) $<

front.o:  front.f90 modules.o
	$(FC) -c   $(FFLAGS) $<

#uzawa.o:  uzawa.f90 modules.o
#	$(FC) -c   $(FFLAGS) $<

%.o : %.f90
	$(FC) -c  $(FFLAGS) $<

compare: compare.o
	@$(CC) -o compare compare.o -lm
	@mv compare ~/bin

parisdeconv: parisdeconv.o
	@$(CC) -o parisdeconv parisdeconv.o -lm
	@mv parisdeconv ~/bin

.c.o:   $< 
	@cc -c $(CFLAGS)   $< 
