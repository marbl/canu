PROG    = test-bitPackedFile \
          test-bitPacking \
          test-freeDiskSpace \
          test-intervalList \
          test-md5 \
          test-mmap \
          test-palloc \
          test-readBuffer \
          test-types

INCLUDE = -I..
LIBS    = -L.. -lutil -lm
OBJS    = 

include ../../Make.compilers

all: $(PROG)
	@echo Tests passed!

test-bitPackedFile: test-bitPackedFile.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bitPackedFile.o test-bitPackedFile.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bitPackedFile test-bitPackedFile.o $(LIBS)
	./test-bitPackedFile 500000 200

test-bitPackedArray: test-bitPackedArray.C ../libutil.a ../util.h ../util++.H
	$(CXX) $(CXXFLAGS_COMPILE) -c -o test-bitPackedArray.o test-bitPackedArray.C $(INCLUDE)
	$(CXX) $(CXXLDFLAGS) -o test-bitPackedArray test-bitPackedArray.o $(LIBS)
	./test-bitPackedArray

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
	./test-readBuffer

test-types: test-types.c ../libutil.a ../util.h ../util++.H
	$(CC) $(CFLAGS_COMPILE) -c -o test-types.o test-types.c $(INCLUDE)
	$(CC) $(CLDFLAGS) -o test-types test-types.o $(LIBS)
	./test-types

clean:
	rm -f $(PROG) *.o *junk*
