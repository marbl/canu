PROG    = 
INCLUDE = -I../libbri
LIBS    = -L../libbri -lbri

OBJS     = positionDB.o positionDB-sort.o positionDB-access.o positionDB-dump.o

all: $(OBJS)

clean:
	rm -f *.o $(PROG)

include ../Makefile-configuration
