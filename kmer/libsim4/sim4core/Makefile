PROG    = sim4dbseq
INCLUDE = -I../libbri
LIBS    = -L../libbri -lbri -lm

SRCS    = Xtend1.C align.C exon.H exonManager.H exon_cores.C extend.C \
          greedy.C log.C msp-link.C msp-sort.C mspManager.H \
          pluri_align.C poly.C sim4.H sim4b1.C sim4db-seq.C sim4db.C sim4db.H \
          sim4dbthr.C sim4defines.H sim4parameters.C sim4parameters.H sim4string.C \
          splice.C table.C util.C
OBJS    = sim4string.o sim4parameters.o Xtend1.o sim4b1.o align.o \
          splice.o greedy.o extend.o util.o exon_cores.o pluri_align.o \
          log.o msp-link.o msp-sort.o table.o poly.o

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE   = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

all:  $(PROG)

sim4dbseq: buildinfo-sim4.o sim4db.o sim4db-seq.o $(OBJS)
	$(CC) $(CCOPTS) -o sim4dbseq sim4db.o sim4db-seq.o buildinfo-sim4.o $(OBJS) $(LIBS) $(EFENCE)

buildinfo-sim4.c: $(SRCS)
	../util/buildinfo.pl sim4 "$(CCOPTS)" $(SRCS)

clean:
	rm -f $(PROG) *.o buildinfo-sim4.*

include ../Makefile-configuration
