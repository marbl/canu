PROG    = snapper2
INCLUDE = -I../existDB -I../positionDB -I../sim4 -I../sim4polish -I../fasta -I../libbri $(THRINC)
LIBS    = -L../existDB -L../positionDB -L../sim4 -L../sim4polish -L../fasta -L../libbri \
          -lexistDB -lpositionDB -lsim4 -lsim4polish -lfasta -lbri

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

SRCS      = snapper2.C configuration.C encodedQuery.C \
            thr-search.C thr-loader.C thr-deadlock.C \
            hitMatrix.C hitMatrix-sort.C \
            hitMatrix.H intervalList.H posix.H snapper2.H
OBJS      = snapper2.o configuration.o encodedQuery.o \
            thr-search.o thr-loader.o thr-deadlock.o hitMatrix-sort.o

all: $(PROG)

#  -static needed for debug on sun with gcc
#  -pthread removed from lines below -- didn't work on solaris

snapper2: hitMatrix.o buildinfo-snapper2.o $(OBJS)
	$(CC) $(CCOPTS) $(INCLUDE) -o snapper2 $(OBJS) hitMatrix.o buildinfo-snapper2.o $(LIBS) $(THRLIB) $(EFENCE) -lm

buildinfo-snapper2.c: $(SRCS)
	../util/buildinfo.pl snapper2 "$(CCOPTS)" $(SRCS)

clean:
	rm -f *.o $(PROG) buildinfo-snapper2.*

include ../Makefile-configuration
