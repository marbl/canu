PROG    = test-bitPackedFile \
          test-bitPacking \
          test-chainedSequence \
          test-intervalList \
          test-md5 \
          test-merstream \
          test-merstreamfile \
          test-merstream-speed \
          test-mmap \
          test-palloc \
          test-readBuffer \
          test-types

INCLUDE = -I..
LIBS    = -L.. -lbri -lm
OBJS    = 

include ../../Make.compilers

all: $(PROG)
	@echo Tests passed!

test-bitPackedFile: test-bitPackedFile.C ../libbri.a ../bri.h ../bri++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bitPackedFile.o test-bitPackedFile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bitPackedFile test-bitPackedFile.o $(LIBS)
	./test-bitPackedFile 500000 200

test-bitPacking: test-bitPacking.C ../libbri.a ../bri.h ../bri++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bitPacking.o test-bitPacking.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bitPacking test-bitPacking.o $(LIBS)
	./test-bitPacking

test-chainedSequence: test-chainedSequence.C ../libbri.a ../bri.h ../bri++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-chainedSequence.o test-chainedSequence.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-chainedSequence test-chainedSequence.o $(LIBS)
	./test-chainedSequence

test-intervalList: test-intervalList.C ../libbri.a ../bri.h ../bri++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-intervalList.o test-intervalList.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-intervalList test-intervalList.o $(LIBS)
	./test-intervalList

test-md5: test-md5.c ../libbri.a ../bri.h ../bri++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-md5.o test-md5.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-md5 test-md5.o $(LIBS)
	./test-md5

test-merstream: test-merstream.C ../libbri.a ../bri.h ../bri++.H ../../fasta/leaff
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream.o test-merstream.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream test-merstream.o $(LIBS)
	./test-merstream
	../../fasta/leaff -G 1000 10000 30000 > junk.fasta
	./test-merstream junk.fasta
	rm -f junk.fasta junk.fastaidx

test-merstreamfile: test-merstreamfile.C ../libbri.a ../bri.h ../bri++.H ../../fasta/leaff
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstreamfile.o test-merstreamfile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstreamfile test-merstreamfile.o $(LIBS)
	../../fasta/leaff -G 1000 5000 10000 > junk.fasta
	./test-merstreamfile junk.fasta
	rm -f junk.fasta junk.fastaidx

test-merstream-speed: test-merstream-speed.C ../libbri.a ../bri.h ../bri++.H ../../fasta/leaff
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-merstream-speed.o test-merstream-speed.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-merstream-speed test-merstream-speed.o $(LIBS)
	../../fasta/leaff -G 10000 1000 10000 > junk.fasta
	cat junk.fasta > /dev/null
	./test-merstream-speed junk.fasta
	rm -f junk.fasta junk.fastaidx

test-mmap: test-mmap.c ../libbri.a ../bri.h ../bri++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-mmap.o test-mmap.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-mmap test-mmap.o $(LIBS)
	./test-mmap 16

test-palloc: test-palloc.c ../libbri.a ../bri.h ../bri++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-palloc.o test-palloc.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-palloc test-palloc.o $(LIBS)
	./test-palloc

test-readBuffer: test-readBuffer.C ../libbri.a ../bri.h ../bri++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-readBuffer.o test-readBuffer.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-readBuffer test-readBuffer.o $(LIBS)
	./test-readBuffer

test-types: test-types.c ../libbri.a ../bri.h ../bri++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-types.o test-types.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-types test-types.o $(LIBS)
	./test-types

clean:
	rm -f $(PROG) *.o *junk*
