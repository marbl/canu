PROG    = sanity
INCLUDE = -I../common
LIBS    = -L../common -lbri -lm

OBJS    = sanity.o

EFENCE  = -L/dev2/adelcher/EFence -lefence
EFENCE  = 

all:  $(PROG)

sanity: $(OBJS)
	$(CC) -o sanity $(OBJS) $(LIBS) $(EFENCE)

test:
	../meryl/meryl -B -s test-seq.fasta -o t -m 20 -t 20
	../meryl/meryl -S -s t -q test-seq.fasta -o tc.out -d -e -c
	../meryl/meryl -S -s t -q test-seq.fasta -o ta.out -d -e -a
	../meryl/meryl -S -s t -q test-seq.fasta -o tp.out -d -e -p

clean:
	rm -f $(PROG) *.o

include ../Makefile-configuration
