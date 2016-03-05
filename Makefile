#--------- paris-ParisSimulator main Makefile --------------------------

# Change the following line to force openmpi to use another compiler
#OMPI_FC=gfortran44

FC = mpif90


FLAGS = -O3 -cpp  -fimplicit-none # -ffpe-summary=invalid,zero,overflow -DOLD_BDRY_COND -fimplicit-none -fbounds-check # -g 
#FLAGS = -O0 -cpp  -fimplicit-none -fimplicit-none -fbounds-check -g 
FLAGSDEBUG = -cpp -g -gstabs -O2  -fimplicit-none # -ffpe-summary=invalid,zero,overflow
HARDFLAGS = -fbounds-check -ffpe-trap=invalid,zero,overflow -Wall # -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal 

CFLAGS = -O # -g -gstabs
BINDIR = $(HOME)/bin

# select option for hypre
# babbage
# HYPRE_DIR = /share/apps/hypre
# Local
SILO_DIR = $(HOME)/cfd/libs/silo-4.10.2

HYPRE_DIR = $(HOME)/cfd/libs/hypre-2.10.0b/src/lib
HYPRE_LIBS =  -L$(HYPRE_DIR) -lHYPRE 
INC_SILO = -I$(SILO_DIR)/include
#DIR_MPI = /usr/local/include
#INC_MPI = -I$(DIR_MPI)
INC = $(INC_SILO) # $(INC_MPI)


ifdef HAVE_VOFI
AFLAGS = -DHAVE_VOFI $(FLAGS) 
VOFI_DIR =  $(HOME)/lib
VOFI_LIBS = -L$(VOFI_DIR) -lvofi
else
AFLAGS = $(FLAGS)
endif

ifdef PHASE_CHANGE
NEWFLAGS = -DPHASE_CHANGE $(AFLAGS)
else
NEWFLAGS = $(AFLAGS)
endif

ifdef HAVE_SILO
FFLAGS = -DHAVE_SILO $(NEWFLAGS) $(INC)
SILO_LIB = -L$(SILO_DIR)/lib -lsilo -lm -lstdc++
else
FFLAGS = $(NEWFLAGS)
endif

LIBS = $(HYPRE_LIBS) $(VOFI_LIBS) $(SILO_LIB)

#------------------------No changes needed beyond this line----------------------------------------------


OBJ = paris.o solids.o modules.o vofmodules.o front.o surface_tension.o lppmodules.o st_testing.o newsolver.o MGsolver.o freesurface.o boiling.o vofnonmodule.o vof_functions.o

SRC = $(wildcard  *.f90) 

install: $(OBJ)
#	@echo compiler is FC = $(FC), mpi override is OMPI_FC = $(OMPI_FC)
	$(FC) -o paris $(FOPTS) $(OBJ) $(LIBS) 
	@if [ ! -d $(BINDIR) ] ; then echo "directory bin does not exist creating it" ; mkdir $(BINDIR) ; fi 
	mv paris $(BINDIR)/paris
	@find .  -name "*.sh" -exec chmod +x  {} \; 

all: tags install pariscompare parisdeconv

clean:
	@rm -fR *.o *.mod paris stats *~ track out* errftc tmp* *.tmp fort.* *.visit core.* *stats*
	@cd Tests; sh ./clean.sh; cd ..
	@cd Documentation; sh ./cleandoc.sh; cd ..
	@find . -type l -exec /bin/rm {} \;

distclean: clean
	@rm -fR  session* *.xml TAGS tags input

test:  install pariscompare parisdeconv pariscompare3D
	@echo "The test suite takes less than 4 minutes on a 4-core intel i7 MacBookPro"
	@cd Tests; chmod +x ./runtests.sh; ./runtests.sh

hardtest:  pariscompare parisdeconv pariscompare3D
	make clean install FLAGS:="$(FLAGSDEBUG) $(HARDFLAGS)"
	@echo "The test suite takes less than 4 minutes on a 4-core intel i7 MacBookPro"
	@cd Tests; chmod +x ./runtests.sh; ./runtests.sh

longtest:  install pariscompare parisdeconv
	@echo "This test takes 12 minutes on a 4-core intel i7 MacBookPro"
	@cd Tests; chmod +x ./runlongtests.sh; ./runlongtests.sh

# single processor test
minitest: install
	@cd Tests/Mini; chmod +x ./run.sh; ./run.sh

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 

paris.o:  paris.f90 solids.o modules.o vofmodules.o front.o surface_tension.o lppmodules.o st_testing.o newsolver.o MGsolver.o freesurface.o boiling.o
	$(FC) -c $(FFLAGS) $<

vofmodules.o: vofmodules.f90 modules.o
	$(FC) -c $(FFLAGS) $<

lppmodules.o: lppmodules.f90 vofmodules.o modules.o
	$(FC) -c $(FFLAGS) $<

surface_tension.o: surface_tension.f90 vofmodules.o modules.o newsolver.o MGsolver.o
	$(FC) -c $(FFLAGS) $<

st_testing.o: st_testing.f90 vofmodules.o modules.o surface_tension.o boiling.o
	$(FC) -c $(FFLAGS) $<

solids.o:  solids.f90 modules.o
	$(FC) -c $(FFLAGS) $<

front.o:  front.f90 modules.o
	$(FC) -c $(FFLAGS) $<

newsolver.o:  newsolver.f90 modules.o MGsolver.o
	$(FC) -c $(FFLAGS) $<

MGsolver.o:  MGsolver.f90 modules.o
	$(FC) -c $(FFLAGS) $<

freesurface.o: freesurface.f90 modules.o
	$(FC) -c $(FFLAGS) $<

boiling.o: boiling.f90 modules.o
	$(FC) -c $(FFLAGS) $<

averages.o: averages.f90 modules.o solids.o vofmodules.o
	$(FC) -c $(FFLAGS) $<

%.o : %.f90
	$(FC) -c $(FFLAGS) $<

pariscompare: compare.o
	@$(CC) -o pariscompare compare.o -lm
	@mv pariscompare $(BINDIR)


pariscompare3D: compare_4cols.o
	$(CC) -o pariscompare3D compare_4cols.o -lm
	mv pariscompare3D $(BINDIR)

parisdeconv: parisdeconv.o
	@$(CC) -o parisdeconv parisdeconv.o -lm
	@mv parisdeconv $(BINDIR)

.c.o:   $< 
	@cc -c $(CFLAGS) $< 
