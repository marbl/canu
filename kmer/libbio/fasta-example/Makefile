PROG    = example
INCLUDE = -I ../fasta        -I../libbri
LIBS    = -L../fasta -lfasta -L../libbri -lbri

SRCS    = example.C
OBJS    = example.o
          

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

all: $(PROG)

example: example.o
	$(CC) $(CCOPTS) $(INCLUDE) -o example example.o $(LIBS) $(EFENCE) -lm

clean:
	rm -f *.o $(PROG) buildinfo*

include ../Makefile-configuration
