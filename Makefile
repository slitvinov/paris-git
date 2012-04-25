#------------------------------------------------
#-------------- Uses openmpi library --------------
#------------------------------------------------

FOPTS   = -O3 
FCOMP   = openmpif90

HYPRE_DIR = /opt/hypre
HYPRE_LIBS =  -L$(HYPRE_DIR)/lib -lHYPRE 


#-----------------------------------------------------------------------------------------------

OBJ     = ftc3d2011.o

run.exe: $(OBJ)
	$(FCOMP) -o run.exe $(FOPTS) $(OBJ) $(HYPRE_LIBS) 

ftc3d2011.o:  ftc3d2011.f90
	$(FCOMP) -c $(FOPTS) ftc3d2011.f90 

clean:
	rm -f *.o *.mod run.exe

