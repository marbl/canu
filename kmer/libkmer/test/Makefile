PROG    = test
INCLUDE = -I.. -I../../libbri
LIBS    = -L.. -L../../libbri -lbri

OBJS     = ../positionDB.o ../positionDB-sort.o ../positionDB-access.o ../positionDB-dump.o

all: $(PROG)

test: test.o $(OBJS)
	$(CC) $(OPTS) $(INCLUDE) -o test test.o $(OBJS) $(LIBS) $(EFENCE)

clean:
	rm -f *.o $(PROG)

include ../../Makefile-configuration
