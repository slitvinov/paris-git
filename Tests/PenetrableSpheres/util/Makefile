SRC =  kozenycarman.c rockread.c 
CFLAGS = -O  # -g -gstabs
BINDIR = $(HOME)/bin

tags:	$(SRC)
# @SZ Create a tags file named TAGS for use by emacs
	@etags $(SRC)
# @SZ Create a tags file named tags for use by vi or textwrangler
# @SZ On MacOS tags and TAGS are identical ! 
# @SZ	ctags paris.f90 


kozenycarman: kozenycarman.o
	@$(CC) -o kozenycarman kozenycarman.o -lm
	mv kozenycarman $(BINDIR)

rockread: rockread.o
	$(CC) -o rockread rockread.o -lm
	mv rockread $(BINDIR)

all: kozenycarman rockread tags

.c.o:   $< 
	cc -c $(CFLAGS)   $< 
