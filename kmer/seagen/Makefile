PROG    = searchGENOME searchGENOMEexactly hitConverter filterEST filterMRNA sortHits
INCLUDE = -I../libbri -I../existDB -I../positionDB $(THRINC)
LIBS    = -L../libbri buildnumber.o ../existDB/existDB.o \
          ../positionDB/positionDB.o ../positionDB/positionDB-sort.o \
          ../positionDB/positionDB-access.o \
          -lbri

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

OBJS      = searchGENOME.o configuration.o encodedQuery.o \
            thr-search.o thr-loader.o thr-deadlock.o hitMatrix.o hitMatrix-sort.o 
OBJSEXACT = searchGENOME.o configuration.o encodedQuery.o \
            thr-search.o thr-loader.o thr-deadlock.o hitMatrix-Exact.o hitMatrix-sort.o 

all: $(PROG);

#  -static needed for debug on sun with gcc
#  -pthread removed from lines below -- didn't work on solaris

searchGENOME: buildnumber $(OBJS)
	$(CC) $(CCOPTS) $(INCLUDE) -o searchGENOME $(OBJS) $(LIBS) $(THRLIB) $(EFENCE)

searchGENOMEexactly: buildnumber $(OBJSEXACT)
	$(CC) $(CCOPTS) $(INCLUDE) -o searchGENOMEexactly $(OBJSEXACT) $(LIBS) $(THRLIB) $(EFENCE)

hitConverter: hitConverter.o
	$(CC) $(CCOPTS) -o hitConverter hitConverter.o
filterByCount: filterByCount.o
	$(CC) $(CCOPTS) -o filterByCount filterByCount.o
filterEST: filterEST.o
	$(CC) $(CCOPTS) -o filterEST filterEST.o
filterMRNA: filterMRNA.o
	$(CC) $(CCOPTS) -o filterMRNA filterMRNA.o
filterByThresh: filterByThresh.o
	$(CC) $(CCOPTS) -o filterByThresh filterByThresh.o
dumpCrapSeqs: dumpCrapSeqs.o encodedQuery.o
	$(CC) $(CCOPTS) -o dumpCrapSeqs dumpCrapSeqs.o encodedQuery.o $(LIBS)
sortHits: sortHits.o
	$(CC) $(CCOPTS) -o sortHits sortHits.o

f: f.o
	$(CC) $(CCOPTS) -o f f.o

clean:
	rm -f *.o $(PROG)

include ../Makefile-configuration
