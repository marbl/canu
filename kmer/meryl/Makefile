PROG    = meryl merle
INCLUDE = -I../libbri
LIBS    = -L../libbri -lbri

EFENCE   = -L/dev2/adelcher/EFence -lefence
EFENCE   = -L/home/walenzbp/projects/EF -lefence
EFENCE   = 

SRCS     = binaryOp.C build.C checkDescr.C dump.C estimate.C \
           mcBucket.C mcBucket.H mcDescription.H merge.C meryl.C meryl.H \
           outputMer.C outputMer.H scan.C unaryOp.C usage.C
OBJS     = meryl.o mcBucket.o outputMer.o estimate.o build.o scan.o merge.o dump.o \
           usage.o checkDescr.o binaryOp.o unaryOp.o

all: $(PROG)

meryl: buildinfo-meryl.o $(OBJS)
	$(CC) $(CCOPTS) $(INCLUDE) -o meryl buildinfo-meryl.o $(OBJS) $(LIBS) $(EFENCE)

merle: buildinfo-merle.o merle.o
	$(CC) $(CCOPTS) $(INCLUDE) -o merle buildinfo-merle.o merle.o $(LIBS) $(EFENCE)

buildinfo-meryl.c: $(SRCS)
	../util/buildinfo.pl meryl "$(CCOPTS)" $(SRCS)

buildinfo-merle.c: merle.C
	../util/buildinfo.pl merle "$(CCOPTS)" merle.C

clean:
	rm -f *.o $(PROG) buildinfo-meryl.* buildinfo-merle.*

include ../Makefile-configuration
