PROG    = alphabet-generate
OBJS    = alphabet.o bitPackedFile.o merstream.o fasta.o util.o fasta-c.o palloc.o sim4reader.o

# AIX needs ar -X 64

bri-util.a: $(OBJS)
	rm -f libbri.a
	ar ruv libbri.a $(OBJS)

#alphabet.o: alphabet.C
#	$(CC) $(CCOPTS) -c -o $*.o $*.C

alphabet.C: alphabet-generate.C
	$(CC) $(CCOPTS) -DMAIN -o alphabet-generate alphabet-generate.C
	./alphabet-generate > alphabet.C
	rm alphabet-generate

clean:
	rm -f libbri.a *.o

include ../Makefile-configuration
