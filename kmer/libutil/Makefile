PROG    = 
INCLUDE = 
LIBS    = 

SRCS    = alphabet.c alphabet.h bit-operations.H bit-packing.H \
          bitPackedFile.C bitPackedFile.H fastastream.H file.h file.c libbri.H \
          libbritypes.h merstream.C merstream.H \
          palloc.c palloc.h readBuffer.H reversecomplement.c \
          splittowords.H time.C time.H \
          memory.c memory.h \
          md5.c md5.h \
          stat.c stat.h
OBJS    = alphabet.o bitPackedFile.o file.o merstream.o palloc.o \
          reversecomplement.o time.o \
          memory.o \
          md5.o ../external/md5lib/md5c.o \
          stat.o

all: libbri.a $(PROG)

install:

libbri.a: buildinfo-libbri.o $(OBJS)
	rm -f libbri.a
	ar $(AROPTS) libbri.a buildinfo-libbri.o $(OBJS)
	$(cc) $(ccOPTS) -DMAIN -o md5 md5.c -L. -lbri
	./md5 && rm md5

alphabet.c: alphabet-generate.c
	$(cc) $(ccOPTS) -DMAIN -o alphabet-generate alphabet-generate.c
	./alphabet-generate
	rm alphabet-generate

buildinfo-libbri.c: $(SRCS)
	../util/buildinfo.pl libbri "$(CCOPTS)" $(SRCS)

clean:
	rm -f libbri.a *.o buildinfo-libbri.*

include ../Makefile-configuration
