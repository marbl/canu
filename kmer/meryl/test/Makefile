PROG    = sanity
INCLUDE = -I.. -I../../libutil
LIBS    = -L.. -L../../libutil -lbio -lutil -lm
OBJS    = sanity.o

include ../../Make.compilers

all:  $(PROG)

sanity: $(OBJS)
	$(CXX) $(CXXFLAGS_COMPILE) -c -o sanity.o sanity.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o sanity sanity.o $(LIBS)
	./sanity

test:
	../meryl -B -s      test-seq1.fasta -o t -m 20


nothin:
	../meryl -S -s t -q test-seq1.fasta -o tc.out -d -e -c
	../meryl -S -s t -q test-seq1.fasta -o ta.out -d -e -a
	../meryl -S -s t -q test-seq1.fasta -o tp.out -d -e -p

clean:
	rm -f $(PROG) *.o
