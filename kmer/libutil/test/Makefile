PROG    = palloc-test reader-delete readBuffer-test
INCLUDE = -I..
LIBS    = -L.. -lbri -lm
OBJS    = 

all: $(PROG)

palloc-test: palloc.c
	$(cc) -DPALLOC_TEST -g -o palloc-test palloc.c

reader-delete: reader-delete.o ../sim4reader.o
	$(CC) -gall -o reader-delete reader-delete.o ../sim4reader.o -lm

readBuffer-test: readBuffer-test.C readBuffer.H
	$(CC) $(CCOPTS) -o readBuffer-test readBuffer-test.C

clean:
	rm -f $(PROG) *.o

include ../../Makefile-configuration
