PROG    = 
INCLUDE = 
LIBS    = 

SRCS    = alphabet.c bit-operations.H bit-packing.H \
          bitPackedFile.C bitPackedFile.H fastastream.H file.h file.c libbri.H \
          libbritypes.h merstream.C merstream.H \
          palloc.c palloc.h readBuffer.H reversecomplement.c \
          sim4reader.c sim4reader.h splittowords.H time.C time.H
OBJS    = alphabet.o bitPackedFile.o file.o merstream.o palloc.o \
          reversecomplement.o sim4reader.o time.o

all: libbri.a $(PROG)

libbri.a: buildinfo-libbri.o $(OBJS)
	rm -f libbri.a
	ar $(AROPTS) libbri.a buildinfo-libbri.o $(OBJS)

alphabet.c: alphabet-generate.c
	$(cc) $(ccOPTS) -DMAIN -o alphabet-generate alphabet-generate.c
	./alphabet-generate
	rm alphabet-generate

buildinfo-libbri.c: $(SRCS)
	../util/buildinfo.pl libbri "$(CCOPTS)" $(SRCS)

clean:
	rm -f libbri.a *.o buildinfo-libbri.*

include ../Makefile-configuration
