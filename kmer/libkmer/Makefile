PROG    = 
INCLUDE = -I../libbri
LIBS    = -L../libbri -lbri

SRCS     = positionDB.C positionDB-access.C positionDB-dump.C positionDB-sort.C positionDB.H
OBJS     = positionDB.o positionDB-access.o positionDB-dump.o positionDB-sort.o

all: buildinfo-positionDB.o $(OBJS)
	rm -f libpositionDB.a
	ar ruv libpositionDB.a buildinfo-positionDB.o $(OBJS)

buildinfo-positionDB.c: $(SRCS)
	../util/buildinfo.pl positionDB "$(CCOPTS)" $(SRCS)

clean:
	rm -f *.o $(PROG) buildinfo-positionDB.*

include ../Makefile-configuration
