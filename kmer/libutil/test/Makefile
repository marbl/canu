PROG    = test-bigQueue \
          test-bitPackedArray \
          test-bitPackedFile \
          test-bitPacking \
          test-freeDiskSpace \
          test-intervalList \
          test-logMsg \
          test-md5 \
          test-mmap \
          test-palloc \
          test-readBuffer \
          test-recordFile \
          test-types \
          tcat

#  Broken, don't test.
#test-bzipBuffer

INCLUDE = -I..
LIBS    = -L.. -lutil -lm
OBJS    = 

include ../../Make.compilers

all: $(PROG)
	@echo Tests passed!

test-bigQueue: test-bigQueue.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bigQueue.o test-bigQueue.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bigQueue test-bigQueue.o $(LIBS)
	time ./test-bigQueue
	time sort -k2n < junk-bigQueue-out-2 > junk-bigQueue-out-2.s
	diff junk-bigQueue-out-1 junk-bigQueue-out-2.s

test-bitPackedArray: test-bitPackedArray.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bitPackedArray.o test-bitPackedArray.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bitPackedArray test-bitPackedArray.o $(LIBS)
	./test-bitPackedArray

test-bitPackedFile: test-bitPackedFile.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bitPackedFile.o test-bitPackedFile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bitPackedFile test-bitPackedFile.o $(LIBS) -lbz2
	./test-bitPackedFile

test-bitPacking: test-bitPacking.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bitPacking.o test-bitPacking.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bitPacking test-bitPacking.o $(LIBS)
	./test-bitPacking

test-freeDiskSpace: test-freeDiskSpace.c ../libutil.a
	$(CC) $(CFLAGS_COMPILE) -c -o test-freeDiskSpace.o test-freeDiskSpace.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-freeDiskSpace test-freeDiskSpace.o $(LIBS)
	./test-freeDiskSpace test-freeDiskSpace.c

test-intervalList: test-intervalList.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-intervalList.o test-intervalList.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-intervalList test-intervalList.o $(LIBS)
	./test-intervalList

test-logMsg: test-logMsg.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-logMsg.o test-logMsg.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-logMsg test-logMsg.o $(LIBS)
	./test-logMsg | wc

test-md5: test-md5.c ../libutil.a ../util.h ../util++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-md5.o test-md5.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-md5 test-md5.o $(LIBS)
	./test-md5

test-mmap: test-mmap.c ../libutil.a ../util.h ../util++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-mmap.o test-mmap.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-mmap test-mmap.o $(LIBS)
	./test-mmap 16

test-palloc: test-palloc.c ../libutil.a ../util.h ../util++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-palloc.o test-palloc.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-palloc test-palloc.o $(LIBS)
	./test-palloc

test-readBuffer: test-readBuffer.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-readBuffer.o test-readBuffer.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-readBuffer test-readBuffer.o $(LIBS)
	./test-readBuffer test-readBuffer

test-recordFile: test-recordFile.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-recordFile.o test-recordFile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-recordFile test-recordFile.o $(LIBS)
	./test-recordFile

test-bzipBuffer: test-bzipBuffer.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bzipBuffer.o test-bzipBuffer.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bzipBuffer test-bzipBuffer.o -lbz2 $(LIBS)
	bzip2 -9vc ./test-bzipBuffer > ./test-bzipBuffer.bz2
	./test-bzipBuffer ./test-bzipBuffer.bz2

test-types: test-types.c ../libutil.a ../util.h ../util++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-types.o test-types.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-types test-types.o $(LIBS)
	./test-types

tcat: tcat.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o tcat.o tcat.C $(INCLUDE) $(THREADS)
	$(CXX) $(CXXLDFLAGS) -o tcat tcat.o $(LIBS) $(THREADL)

clean:
	rm -f $(PROG) *.o *junk*
