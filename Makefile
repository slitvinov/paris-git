#------------------------------------------------
#-------------- Uses openmpi library --------------
#------------------------------------------------

FOPTS   = -O3 
FCOMP   = mpif90

HYPRE_DIR = /home/zaleski/hypre-2.8.0b/src
HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 


#-----------------------------------------------------------------------------------------------

OBJ     = ftc3d2011.o

run.exe: $(OBJ)
	$(FCOMP) -o run.exe $(FOPTS) $(OBJ) $(HYPRE_LIBS) 
#	put it in a sound location @SZ
	mv run.exe ~/bin/ftc3d2011


ftc3d2011.o:  ftc3d2011.f90
	$(FCOMP) -c $(FOPTS) ftc3d2011.f90 
# @SZ Create a tags file named TAGS for use by emacs
	etags ftc3d2011.f90 
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags ftc3d2011.f90 

clean:
	rm -fR *.o *.mod run.exe *.gz stats *~ track out
	cd Speed_Measurement; make clean; cd ..
	cd Poiseuille_Test; make clean; cd ..


