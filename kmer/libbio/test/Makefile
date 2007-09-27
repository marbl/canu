

#PROG    = test-merstream-from-seqstore

PROG    = dump-merstreamfile \
          test-seqStream \
          test-chainedSequence \
          test-fasta-accessor \
          test-merstream \
          test-merstreamfile \
          test-setbits \
          halign-test

DEAD =    test-merstream-speed \
          test-bigmer-msf \


INCLUDE = -I.. -I../../libutil
LIBS    = -L.. -L../../libutil -lbio -lutil -lm
OBJS    = 

include ../../Make.compilers

all: $(PROG)
	@echo Tests passed!

dump-merstreamfile: dump-merstreamfile.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o dump-merstreamfile.o dump-merstreamfile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o dump-merstreamfile dump-merstreamfile.o $(LIBS)

test-merstream-from-seqstore: test-merstream-from-seqstore.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream-from-seqstore.o test-merstream-from-seqstore.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream-from-seqstore test-merstream-from-seqstore.o $(LIBS)

test-seqStream: test-seqStream.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-seqStream.o test-seqStream.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-seqStream test-seqStream.o $(LIBS)
	../../leaff/leaff -G 3 30 40 > junk2.fasta
	./test-seqStream junk2.fasta
	rm -f junk*

test-chainedSequence: test-chainedSequence.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-chainedSequence.o test-chainedSequence.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-chainedSequence test-chainedSequence.o $(LIBS)
	./test-chainedSequence
	../../leaff/leaff -G 1000 1000 3000 > junk2.fasta
	./test-chainedSequence junk2.fasta
	rm -f junk*

test-fasta-accessor: test-fasta-accessor.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-fasta-accessor.o test-fasta-accessor.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-fasta-accessor test-fasta-accessor.o $(LIBS)
	./test-fasta-accessor

test-merstream: test-merstream.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream.o test-merstream.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream test-merstream.o $(LIBS)
	../../leaff/leaff -G 1 10000 30000 > junk.fasta
	./test-merstream junk.fasta
	rm -f junk*
	../../leaff/leaff -G 1000 10000 30000 > junk.fasta
	./test-merstream junk.fasta
	rm -f junk*

test-merstreamfile: test-merstreamfile.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstreamfile.o test-merstreamfile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstreamfile test-merstreamfile.o $(LIBS)
	#../../leaff/leaff -G 2 50 50 > junk.fasta
	../../leaff/leaff -G 100000 10 600 > junk.fasta
	#../../leaff/leaff -G 10000 10 10000 > junk.fasta
	#../../leaff/leaff -G 30000 10000 10000 > junk.fasta
	./test-merstreamfile junk.fasta
	#rm -f junk.fasta junk.fastaidx

test-merstream-speed: test-merstream-speed.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream-speed.o test-merstream-speed.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream-speed test-merstream-speed.o $(LIBS)
	../../leaff/leaff -G 10000 1000 10000 > junk.fasta
	cat junk.fasta > /dev/null
	./test-merstream-speed junk.fasta
	rm -f junk*

test-setbits: test-setbits.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-setbits.o test-setbits.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-setbits test-setbits.o $(LIBS)
	./test-setbits

test-bigmer-msf: test-bigmer-msf.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bigmer-msf.o test-bigmer-msf.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bigmer-msf test-bigmer-msf.o $(LIBS)
	./test-bigmer-msf

halign-test: halign-test.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o halign-test.o halign-test.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o halign-test halign-test.o $(LIBS)

clean:
	rm -f $(PROG) *.o *junk*
