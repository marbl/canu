PROG    = filterPolishes mergePolishes stripPolishes sortPolishes \
          pickBestPolish pickBestPolish-debug parseSNP \
          plotCoverageVsIdentity cleanPolishes uniqPolishes
INCLUDE = -I../libbri
LIBS    = -L../libbri -lbri -lm

OBJS    =

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

all: $(PROG)


filterPolishes: filterPolishes.o $(OBJS)
	$(CC) $(CCOPTS) -o filterPolishes filterPolishes.o $(OBJS) $(LIBS) $(EFENCE)

mergePolishes: mergePolishes.o $(OBJS)
	$(CC) $(CCOPTS) -o mergePolishes mergePolishes.o $(OBJS) $(LIBS) $(EFENCE)

stripPolishes: stripPolishes.o $(OBJS)
	$(CC) $(CCOPTS) -o stripPolishes stripPolishes.o $(OBJS) $(LIBS) $(EFENCE)

sortPolishes: sortPolishes.o $(OBJS)
	$(CC) $(CCOPTS) -o sortPolishes sortPolishes.o $(OBJS) $(LIBS) $(EFENCE)

pickBestPolish: pickBestPolish.o $(OBJS)
	$(CC) $(CCOPTS) -o pickBestPolish pickBestPolish.o $(OBJS) $(LIBS) $(EFENCE)
pickBestPolish-debug: pickBestPolish-debug.o $(OBJS)
	$(CC) $(CCOPTS) -o pickBestPolish-debug pickBestPolish-debug.o $(OBJS) $(LIBS) $(EFENCE)

pickBestPolish-debug.o: pickBestPolish.c
	$(cc) $(ccOPTS) $(WARNS) $(INCLUDE) -DVALIDATE -c -o pickBestPolish-debug.o pickBestPolish.c

cleanPolishes: cleanPolishes.o $(OBJS)
	$(CC) $(CCOPTS) -o cleanPolishes cleanPolishes.o $(OBJS) $(LIBS) $(EFENCE)
plotIntronSize: plotIntronSize.o $(OBJS)
	$(CC) $(CCOPTS) -o plotIntronSize plotIntronSize.o $(OBJS) $(LIBS) $(EFENCE)

plotCoverageVsIdentity: plotCoverageVsIdentity.o $(OBJS)
	$(CC) $(CCOPTS) -o plotCoverageVsIdentity plotCoverageVsIdentity.o $(OBJS) $(LIBS) $(EFENCE)

parseSNP: parseSNP.o $(OBJS)
	$(CC) $(CCOPTS) -o parseSNP parseSNP.o $(OBJS) $(LIBS) $(EFENCE)

#  Compare two polished sets
comparePolishes: comparePolishes.o $(OBJS)
	$(CC) $(CCOPTS) -o comparePolishes comparePolishes.o $(OBJS) $(LIBS) $(EFENCE)

#  Trim sequences based on the largest region mapped
trimSequencesBasedOnMatches: trimSequencesBasedOnMatches.o $(OBJS)
	$(CC) $(CCOPTS) -o trimSequencesBasedOnMatches trimSequencesBasedOnMatches.o $(OBJS) $(LIBS) $(EFENCE)

uniqPolishes: uniqPolishes.o $(OBJS)
	$(CC) $(CCOPTS) -o uniqPolishes uniqPolishes.o $(OBJS) $(LIBS) $(EFENCE)

clean:
	rm -f $(PROG) *.o

include ../Makefile-configuration
