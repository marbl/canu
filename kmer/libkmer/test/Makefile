include ../../Make.compilers

passed: tester ../../fasta/leaff
	../../fasta/leaff -G 1000 5000 10000 > x.fasta
	./tester -test1 x.fasta
	../../fasta/leaff -G 100 500 1000 > y.fasta
	./tester -test2 x.fasta y.fasta
	rm -f x.fasta y.fasta
	touch passed


INCLUDE = -I../../libbri -I../../positionDB               -I../../existDB 
LIBS    = -L../../libbri -L../../positionDB -L../../meryl -L../../existDB \
          -lpositionDB -lexistDB -lmeryl -lbri
DBGOPT  = -DERROR_CHECK_COUNTING -DERROR_CHECK_COUNTING_ENCODING -DERROR_CHECK_EMPTY_BUCKETS


tester: ../positionDB.C ../positionDB.H ../positionDB-access.C ../positionDB-dump.C ../positionDB-sort.C ../positionDB-file.C
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o driver.o            ../driver.C            $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB.o        ../positionDB.C        $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-access.o ../positionDB-access.C $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-dump.o   ../positionDB-dump.C   $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-sort.o   ../positionDB-sort.C   $(INCLUDE)
	$(CXX) $(DBGOPT) $(CXXFLAGS_COMPILE) -c -o positionDB-file.o   ../positionDB-file.C   $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o tester driver.o positionDB.o positionDB-access.o positionDB-dump.o positionDB-sort.o positionDB-file.o $(LIBS)



clean:
	rm -f passed tester *.o x.fasta y.fasta
