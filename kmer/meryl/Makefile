PROG    = meryl
INCLUDE = -I../libbri
LIBS    = -L../libbri -lbri buildnumber.o

EFENCE   = -L/dev2/adelcher/EFence -lefence
EFENCE   = -L/home/walenzbp/projects/EF -lefence
EFENCE   = 

OBJS     = mcBucket.o outputMer.o estimate.o build.o scan.o merge.o dump.o usage.o \
           checkDescr.o binaryOp.o unaryOp.o

all: meryl

meryl: buildnumber meryl.o $(OBJS)
	$(CC) $(CCOPTS) $(INCLUDE) -o meryl meryl.o $(OBJS) $(LIBS) $(EFENCE)

clean:
	rm -f *.o $(PROG)

include ../Makefile-configuration
