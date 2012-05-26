#------------------------------------------------
#-------------- Uses openmpi library --------------
#------------------------------------------------

FFLAGS   = -O3 
FC   = mpif90
CFLAGS = -O

#HYPRE_DIR = /home/zaleski/hypre-2.8.0b/src
HYPRE_DIR = /opt/hypre
HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 


#-----------------------------------------------------------------------------------------------

OBJ     = ftc3d2011.o utilc.o solids.o modules1.o
SRC     = ftc3d2011.f90 utilc.c solids.f90 modules1.f90

install: $(OBJ)
	$(FC) -o run.exe $(FOPTS) $(OBJ) $(HYPRE_LIBS) 
#	put it in a sound location @SZ
	mv run.exe ~/bin/ftc3d2011


ftc3d2011.o:  ftc3d2011.f90 solids.o modules1.o
	$(FC) -c $(FFLAGS) ftc3d2011.f90 

modules: solids.o modules1.o

all: tags install

clean:
	rm -fR *.o *.mod run.exe *.gz stats *~ track out
	cd Speed_Measurement; make clean; cd ..
	cd Poiseuille_Test; make clean; cd ..

solids.o:  solids.f90 modules1.o
	$(FC) -c $(FFLAGS) $<


%.o : %.f90
	$(FC) -c $(FFLAGS) $<

.c.o:   $< utilc.h
	$(CC) -c $(CFLAGS)   $< 

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags ftc3d2011.f90 

