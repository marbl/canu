PROG    = alphabet-generate
INCLUDE = 
LIBS    = 

SRCS    = alphabet-generate.C alphabet.C bit-packing.H bitPackedFile.C bitPackedFile.H \
          fasta-c.C fasta.C fasta.H fastastream.H libbri.H libbritypes.H \
          merstream.C merstream.H outputString.H sequencemanager-compressed.H \
          sequencemanager.H speedCounter.H util.C
OBJS    = alphabet.o bitPackedFile.o merstream.o fasta.o util.o fasta-c.o palloc.o sim4reader.o

# AIX needs ar -X 64

libbri.a: buildinfo-libbri.o $(OBJS)
	rm -f libbri.a
	ar ruv libbri.a buildinfo-libbri.o $(OBJS)

#alphabet.o: alphabet.C
#	$(CC) $(CCOPTS) -c -o $*.o $*.C

alphabet.C: alphabet-generate.C
	$(CC) $(CCOPTS) -DMAIN -o alphabet-generate alphabet-generate.C
	./alphabet-generate > alphabet.C
	rm alphabet-generate

buildinfo-libbri.c: $(SRCS)
	../util/buildinfo.pl libbri "$(CCOPTS)" $(SRCS)

clean:
	rm -f libbri.a *.o buildinfo-libbri.*

include ../Makefile-configuration
