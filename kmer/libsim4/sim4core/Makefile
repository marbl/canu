PROG    = 
INCLUDE = -I../fasta         -I../sim4polish              -I../libbri
LIBS    = -L../fasta -lfasta -L../sim4polish -lsim4polish -L../libbri -lbri -lm

SRCS    = sim4command.C \
          sim4parameters.C \
          sim4string.C \
          Xtend1.C align.C exon.H exonManager.H exon_cores.C extend.C \
          greedy.C log.C msp-link.C msp-sort.C mspManager.H \
          pluri_align.C poly.C sim4.H \
          sim4b1.C \
          sim4b1a.C \
          sim4b1-1.C \
          sim4b1-2.C \
          sim4b1-3.C \
          sim4b1-4.C \
          sim4defines.H sim4parameters.C sim4parameters.H sim4string.C \
          splice.C table.C util.C
OBJS    = sim4command.o \
          sim4parameters.o \
          sim4string.o \
          Xtend1.o \
          sim4b1.o \
          sim4b1a.o \
          sim4b1-1.o \
          sim4b1-2.o \
          sim4b1-3.o \
          sim4b1-4.o \
          align.o \
          splice.o greedy.o extend.o util.o exon_cores.o pluri_align.o \
          log.o msp-link.o msp-sort.o table.o poly.o

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE   = -L/home/walenzbp/projects/EF -lefence
EFENCE  = 

all: libsim4.a

install:

libsim4.a: $(OBJS) buildinfo-libsim4.o
	rm -f libsim4.a
	ar $(AROPTS) libsim4.a buildinfo-libsim4.o $(OBJS)

buildinfo-libsim4.c: $(SRCS)
	../util/buildinfo.pl libsim4 "$(CCOPTS)" $(SRCS)

clean:
	rm -f $(PROG) *.o *.a buildinfo-*

include ../Makefile-configuration
