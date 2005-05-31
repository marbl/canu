PROG    = stupidcount exhaustive
INCLUDE = -I.. -I../../libutil -I../../libbio -I../../libmeryl
LIBS    = -L.. -L../../libutil -L../../libbio -L../../libmeryl -lmeryl -lbio -lutil -lm

include ../../Make.compilers

all:  $(PROG)

stupidcount: stupidcount.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o stupidcount.o stupidcount.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o stupidcount stupidcount.o $(LIBS)

exhaustive: exhaustive.C kmerlite.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o exhaustive.o exhaustive.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o exhaustive exhaustive.o $(LIBS)

test-exhaustive: exhaustive ../meryl ../../leaff/leaff
	../../leaff/leaff -G 1000 10000 40000 > /scratch/g.fasta
	../meryl -B -s /scratch/g.fasta -o /scratch/s -m 100
	./exhaustive -m /scratch/s -f /scratch/g.fasta

test-reduce: ../meryl
	#  Build the initial table
	../meryl -B -f -m 20 -s test-seq1.fasta -o 1
	#  Dump the initial table as fasta
	../meryl -Dt -n 0 -s 1 > 2.reduce.fasta
	#  Build a new table on the dumped fasta
	../meryl -B -f -m 20 -s 2.reduce.fasta -o 2
	#  Remove one copy of each mer
	../meryl -M sub -s 1 -s 2 -o 3
	#  Dump the resulting file
	../meryl -Dt -n 1 -s 3
	touch test-reduce

test:
	../meryl -B -s      test-seq1.fasta -o t -m 20

nothin:
	../meryl -S -s t -q test-seq1.fasta -o tc.out -d -e -c
	../meryl -S -s t -q test-seq1.fasta -o ta.out -d -e -a
	../meryl -S -s t -q test-seq1.fasta -o tp.out -d -e -p

clean:
	rm -f $(PROG) *.o
