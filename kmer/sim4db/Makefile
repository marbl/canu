PROG    = sim4db
INCLUDE = -I../fasta         -I../sim4        -I../sim4polish              -I../libbri
LIBS    = -L../fasta -lfasta -L../sim4 -lsim4 -L../sim4polish -lsim4polish -L../libbri -lbri -lm

SRCS    = sim4db.C
OBJS    = sim4db.o

EFENCE  = -L/home/walenzbp/projects/EF -lefence
EFENCE  = -L/usr/local/lib -lefence
EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  =

all:  $(PROG)

install: $(PROG)
	cp -p $(PROG) ../bin

#sim4dbg: buildinfo-sim4db.o $(OBJS)
#	$(CC) $(CCOPTS) -o sim4dbg buildinfo-sim4db.o $(OBJS) $(LIBS) $(EFENCE)

sim4db: buildinfo-sim4db.o $(OBJS)
	$(CC) $(CCOPTS) -o sim4db buildinfo-sim4db.o $(OBJS) $(LIBS) $(EFENCE)

buildinfo-sim4db.c: $(SRCS)
	../util/buildinfo.pl sim4db "$(CCOPTS)" $(SRCS)

clean:
	rm -f $(PROG) *.o buildinfo-*

include ../Makefile-configuration
