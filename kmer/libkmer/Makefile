PROG    = 
INCLUDE = -I../libbri
LIBS    = 

SRCS     = positionDB.C positionDB-access.C positionDB-dump.C positionDB-sort.C positionDB.H
OBJS     = positionDB.o positionDB-access.o positionDB-dump.o positionDB-sort.o

all: libpositionDB.a

install:

libpositionDB.a: buildinfo-positionDB.o $(OBJS)
	rm -f libpositionDB.a
	ar $(AROPTS) libpositionDB.a buildinfo-positionDB.o $(OBJS)

buildinfo-positionDB.c: $(SRCS)
	../util/buildinfo.pl positionDB "$(CCOPTS)" $(SRCS)

clean:
	rm -f *.o $(PROG) buildinfo-positionDB.* libpositionDB.a

include ../Makefile-configuration
