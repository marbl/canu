PROG    = reader-delete
INCLUDE = -I..
LIBS    = -L.. -lbri -lm
OBJS    = 

all: $(PROG)

fasta-test: fasta-test.o
	$(CC) $(OPTS) -o fasta-test fasta-test.o $(LIBS)

intervalList-test: intervalList-test.o
	$(CC) $(OPTS) -o intervalList-test intervalList-test.o $(LIBS)

palloc-test: palloc.c
	$(cc) -DPALLOC_TEST -g -o palloc-test palloc.c

reader-delete: reader-delete.o ../sim4reader.o
	$(CC) -gall -o reader-delete reader-delete.o ../sim4reader.o -lm

clean:
	rm -f $(PROG) *.o

include ../../Makefile-configuration
