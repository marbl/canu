PROG    = test-chainedSequence \
          test-fasta-accessor \
          test-merstream \
          test-merstreamfile \
          test-merstream-speed

INCLUDE = -I.. -I../../libutil
LIBS    = -L.. -L../../libutil -lbio -lutil -lm
OBJS    = 

include ../../Make.compilers

all: $(PROG)
	@echo Tests passed!

test-chainedSequence: test-chainedSequence.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-chainedSequence.o test-chainedSequence.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-chainedSequence test-chainedSequence.o $(LIBS)
	./test-chainedSequence

test-fasta-accessor: test-fasta-accessor.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-fasta-accessor.o test-fasta-accessor.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-fasta-accessor test-fasta-accessor.o $(LIBS)
	./test-fasta-accessor

test-merstream: test-merstream.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream.o test-merstream.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream test-merstream.o $(LIBS)
	./test-merstream
	../../leaff/leaff -G 1000 10000 30000 > junk.fasta
	./test-merstream junk.fasta
	rm -f junk.fasta junk.fastaidx

test-merstreamfile: test-merstreamfile.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstreamfile.o test-merstreamfile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstreamfile test-merstreamfile.o $(LIBS)
	../../leaff/leaff -G 1000 5000 10000 > junk.fasta
	./test-merstreamfile junk.fasta
	rm -f junk.fasta junk.fastaidx

test-merstream-speed: test-merstream-speed.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream-speed.o test-merstream-speed.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream-speed test-merstream-speed.o $(LIBS)
	../../leaff/leaff -G 10000 1000 10000 > junk.fasta
	cat junk.fasta > /dev/null
	./test-merstream-speed junk.fasta
	rm -f junk.fasta junk.fastaidx

clean:
	rm -f $(PROG) *.o *junk*
