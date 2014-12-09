

PROG    = test-merstream-speed

INCLUDE = -I.. -I../../libutil -I../../libbio -I../../libseq
LIBS    = -L.. -L../../libutil -L../../libbio -L../../libseq -lseq -lbio -lutil -lm
OBJS    = 

include ../../Make.compilers

all: $(PROG)
	@echo Tests passed!

test-merstream-speed: test-merstream-speed.C
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream-speed.o test-merstream-speed.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream-speed test-merstream-speed.o $(LIBS)
	../../leaff/leaff -G 10000 1000 10000 > junk.fasta
	cat junk.fasta > /dev/null
	./test-merstream-speed junk.fasta
	rm -f junk*

clean:
	rm -f $(PROG) *.o *junk*
