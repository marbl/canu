PROG    = searchGENOME searchGENOMEexactly hitConverter filterEST filterMRNA sortHits
INCLUDE = -I../existDB -I../positionDB -I../fasta -I../libbri $(THRINC)
LIBS    = -L../existDB -L../positionDB -L../fasta -L../libbri \
          -lexistDB -lpositionDB -lfasta -lbri

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

SRCS      = searchGENOME.C configuration.C encodedQuery.C \
            thr-search.C thr-loader.C thr-deadlock.C \
            hitMatrix.C hitMatrix-Exact.C hitMatrix-sort.C \
            aHit.H hitMatrix.H intervalList.H posix.H searchGENOME.H
OBJS      = searchGENOME.o configuration.o encodedQuery.o \
            thr-search.o thr-loader.o thr-deadlock.o hitMatrix-sort.o

all: $(PROG)

#  -static needed for debug on sun with gcc
#  -pthread removed from lines below -- didn't work on solaris

searchGENOME: hitMatrix.o buildinfo-searchGENOME.o $(OBJS)
	$(CC) $(CCOPTS) $(INCLUDE) -o searchGENOME $(OBJS) hitMatrix.o buildinfo-searchGENOME.o $(LIBS) $(THRLIB) $(EFENCE) -lm

searchGENOMEexactly: hitMatrix-Exact.o buildinfo-searchGENOME.o $(OBJS)
	$(CC) $(CCOPTS) $(INCLUDE) -o searchGENOMEexactly $(OBJS) hitMatrix-Exact.o buildinfo-searchGENOME.o $(LIBS) $(THRLIB) $(EFENCE) -lm

hitConverter: hitConverter.o
	$(CC) $(CCOPTS) -o hitConverter hitConverter.o
filterEST: filterEST.o
	$(CC) $(CCOPTS) -o filterEST filterEST.o
filterMRNA: filterMRNA.o
	$(CC) $(CCOPTS) -o filterMRNA filterMRNA.o
sortHits: sortHits.o
	$(CC) $(CCOPTS) -o sortHits sortHits.o

filtertest: filtertest.o
	$(CC) $(CCOPTS) -o filtertest filtertest.o

buildinfo-searchGENOME.c: $(SRCS)
	../util/buildinfo.pl searchGENOME "$(CCOPTS)" $(SRCS)

clean:
	rm -f *.o $(PROG) buildinfo-searchGENOME.*

include ../Makefile-configuration
