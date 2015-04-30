SRC =  kozenycarman.c rockread.c
CFLAGS = -O  # -g -gstabs
BINDIR = $(HOME)/bin
FC = mpif90
FLAGS = -O2 -cpp -fimplicit-none #-fbounds-check # -g 

SILO_DIR = $(HOME)/cfd/libs/silo-4.10.2
INC_SILO = -I$(SILO_DIR)/include
INC = $(INC_SILO)

ifdef HAVE_SILO
FFLAGS = -DHAVE_SILO $(AFLAGS) $(INC)
SILO_LIB = -L$(SILO_DIR)/lib -lsilo -lm -lstdc++
endif

LIBS = $(HYPRE_LIBS) $(VOFI_LIBS) $(SILO_LIB)

all: kozenycarman rockread tags

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 

kozenycarman: kozenycarman.o
	$(CC) -o kozenycarman kozenycarman.o -lm
	mv kozenycarman $(BINDIR)

rockread: rockread.o
	$(CC) -o rockread rockread.o -lm
	

 
post_utility: post_utility.o
	$(FC) -o post_utility post_utility.o $(LIBS)
	mv post_utility $(BINDIR)

post_utility.o: post_utility.f90 
	$(FC) -c $(INC) $<
clean:
	rm *~ *.o

.c.o:   $< 
	cc -c $(CFLAGS)   $< 
	
%.o : %.f90
	$(FC) -c $(FFLAGS) $<
