PROG    = searchGENOME searchGENOMEexactly hitConverter filterEST filterMRNA sortHits
INCLUDE = -I../libbri -I../existDB -I../positionDB $(THRINC)
LIBS    = -L../libbri -L../existDB -L../positionDB \
          -lexistDB -lpositionDB -lbri

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

SRCS      = searchGENOME.C configuration.C encodedQuery.C \
            thr-search.C thr-loader.C thr-deadlock.C \
            hitMatrix.C hitMatrix-Exact.C hitMatrix-sort.C \
            aHit.H hitMatrix.H intervalList.H posix.H searchGENOME.H
OBJS      = searchGENOME.o configuration.o encodedQuery.o \
            thr-search.o thr-loader.o thr-deadlock.o hitMatrix-sort.o

all: buildinfo-searchGENOME.c $(PROG)

#  -static needed for debug on sun with gcc
#  -pthread removed from lines below -- didn't work on solaris

searchGENOME: buildinfo-searchGENOME.o $(OBJS) hitMatrix.o
	$(CC) $(CCOPTS) $(INCLUDE) -o searchGENOME $(OBJS) hitMatrix.o buildinfo-searchGENOME.o $(LIBS) $(THRLIB) $(EFENCE)

searchGENOMEexactly: buildinfo-searchGENOME.o $(OBJS) hitMatrix-Exact.o
	$(CC) $(CCOPTS) $(INCLUDE) -o searchGENOMEexactly $(OBJS) hitMatrix-Exact.o buildinfo-searchGENOME.o $(LIBS) $(THRLIB) $(EFENCE)

hitConverter: hitConverter.o
	$(CC) $(CCOPTS) -o hitConverter hitConverter.o
filterEST: filterEST.o
	$(CC) $(CCOPTS) -o filterEST filterEST.o
filterMRNA: filterMRNA.o
	$(CC) $(CCOPTS) -o filterMRNA filterMRNA.o
sortHits: sortHits.o
	$(CC) $(CCOPTS) -o sortHits sortHits.o

buildinfo-searchGENOME.c: $(SRCS)
	../util/buildinfo.pl searchGENOME "$(CCOPTS)" $(SRCS)

clean:
	rm -f *.o $(PROG) buildinfo-searchGENOME.*

include ../Makefile-configuration
