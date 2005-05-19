include ../../Make.compilers

#  Bigger tblsize makes existDB much faster, but uses more memory (not
#  much, really).  23 is nice.

all: position-passed
	@echo "existDB has expensive tests.  They take:"
	@echo "  70 minutes on 2.0GHz G5 (old)"
	@echo "  60 minutes on 2.8GHz P4 (with KMER=1)"
	@echo " 120 minutes on 2.0GHz G5 (with KMER=8)"
	@echo "If you really want to run them, do 'make exist-passed'."

position-passed: position-passed1 position-passed2 exist-fast
	touch position-passed
xp.fasta: ../../leaff/leaff
	../../leaff/leaff -G 1000 5000 10000 > xp.fasta
yp.fasta: ../../leaff/leaff
	../../leaff/leaff -G 100 500 1000 > yp.fasta
position-passed1: testerp xp.fasta
	./testerp -test1 xp.fasta
	touch position-passed1
position-passed2: testerp xp.fasta yp.fasta
	./testerp -test2 xp.fasta yp.fasta
	touch position-passed2
exist-fast-passed: ../existDB xe.fasta
	../existDB -mersize 14 -tblsize 21 -testfiles xe.fasta junk
	../existDB -mersize 14 -tblsize 21 -testexistence xe.fasta
	../../meryl/meryl -B -f -m 14 -s xe.fasta -o xe
	../existDB -mersize 14 -tblsize 21 -testexhaustive xe.fasta xe
	touch exist-fast-passed


exist-passed: exist-passed1 exist-passed2 exist-passed3
	touch exist-passed
xe.fasta: ../../leaff/leaff
	../../leaff/leaff -G 1000 5000 10000 > xe.fasta
exist-passed1: ../existDB xe.fasta
	../existDB -mersize 17 -tblsize 23 -testfiles xe.fasta junk
	touch exist-passed1
exist-passed2: ../existDB xe.fasta
	../existDB -mersize 17 -tblsize 23 -testexistence xe.fasta
	touch exist-passed2
exist-passed3: ../existDB xe.fasta ../../meryl/meryl
	../../meryl/meryl -B -f -m 17 -s xe.fasta -o xe
	../existDB -mersize 17 -tblsize 23 -testexhaustive xe.fasta xe
	touch exist-passed3


INCLUDE = -I../../libbio -I../../libutil -I../../libmeryl
LIBS    = -L../../libbio -L../../libutil -L../../libmeryl -L.. -lkmer -lmeryl -lbio -lutil
DBGOPT  = -DERROR_CHECK_COUNTING -DERROR_CHECK_COUNTING_ENCODING -DERROR_CHECK_EMPTY_BUCKETS

testerp: ../positionDB.C ../positionDB.H ../positionDB-access.C ../positionDB-dump.C ../positionDB-sort.C ../positionDB-file.C
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o driverp.o           ../driver-posDB.C      $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB.o        ../positionDB.C        $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-access.o ../positionDB-access.C $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-dump.o   ../positionDB-dump.C   $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-sort.o   ../positionDB-sort.C   $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-file.o   ../positionDB-file.C   $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o testerp driverp.o positionDB.o positionDB-access.o positionDB-dump.o positionDB-sort.o positionDB-file.o $(LIBS) -lm

#  XXX:  There isn't any reason we need to build testere, we could
#  just use ../existDB (as it did before!)
testere: ../existDB.C ../existDB-create-from-fasta.C ../existDB-create-from-meryl.C ../existDB-state.C
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o drivere.o                    ../driver-existDB.C             $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o existDB.o                    ../existDB.C                    $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o existDB-create-from-fasta.o  ../existDB-create-from-fasta.C  $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o existDB-create-from-meryl.o  ../existDB-create-from-meryl.C  $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o existDB-state.o              ../existDB-state.C              $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o testere drivere.o existDB.o existDB-create-from-fasta.o existDB-create-from-meryl.o existDB-state.o $(LIBS) -lm

clean:
	rm -f *passed*
	rm -f testerp *.o xp.fasta* yp.fasta*
	rm -f testere junk* xe.fasta* xe.mcidx xe.mcdat xe.merStream
