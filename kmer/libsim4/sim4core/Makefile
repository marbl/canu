PROG    = findcdnagapped sim4dbseq sim4dbthr
INCLUDE = -I../libbri
LIBS    = -L../libbri buildnumber.o -lbri -lm

OBJS    = sim4string.o sim4parameters.o Xtend1.o sim4b1.o align.o \
          splice.o greedy.o extend.o util.o exon_cores.o pluri_align.o \
          log.o msp-link.o msp-sort.o table.o poly.o

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = 

#all:  $(PROG)
all:  sim4dbseq findcdnagapped

#  For alpha, remove -pthread and add -lpthread -lrt to the link line.
#sim4dbthr: buildnumber sim4db.o sim4db-thr.o $(OBJS)
#	$(CC) $(CCOPTS) -o sim4dbthr sim4db.o sim4db-thr.o $(OBJS) $(LIBS) -lpthread -lrt $(EFENCE)

sim4dbseq: buildnumber sim4db.o sim4db-seq.o $(OBJS)
	$(CC) $(CCOPTS) -o sim4dbseq sim4db.o sim4db-seq.o $(OBJS) $(LIBS) $(EFENCE)

findcdnagapped: buildnumber findcdnagapped.o $(OBJS)
	$(CC) $(CCOPTS) -o findcdnagapped findcdnagapped.o $(OBJS) $(LIBS) $(EFENCE)

clean:
	rm -f $(PROG) *.o

include ../Makefile-configuration
