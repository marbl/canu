PROG    = readBuffer-test
INCLUDE = 
LIBS    = 

SRCS    = alphabet-generate.C alphabet.C bit-packing.H bitPackedFile.C bitPackedFile.H \
          fastastream.H libbri.H libbritypes.H \
          merstream.C merstream.H outputString.H \
          speedCounter.H util.C readBuffer.H readBuffer-test.C
OBJS    = alphabet.o bitPackedFile.o merstream.o util.o palloc.o sim4reader.o

all: libbri.a $(PROG)

libbri.a: buildinfo-libbri.o $(OBJS)
	rm -f libbri.a
	ar $(AROPTS) libbri.a buildinfo-libbri.o $(OBJS)

#alphabet.o: alphabet.C
#	$(CC) $(CCOPTS) -c -o $*.o $*.C

alphabet.C: alphabet-generate.C
	$(CC) $(CCOPTS) -DMAIN -o alphabet-generate alphabet-generate.C
	./alphabet-generate > alphabet.C
	rm alphabet-generate

readBuffer-test: readBuffer-test.C readBuffer.H
	$(CC) $(CCOPTS) -o readBuffer-test readBuffer-test.C

buildinfo-libbri.c: $(SRCS)
	../util/buildinfo.pl libbri "$(CCOPTS)" $(SRCS)

clean:
	rm -f libbri.a *.o buildinfo-libbri.*

include ../Makefile-configuration
