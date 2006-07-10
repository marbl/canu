
# Ross's build system wasn't doing well here.  We need to include
# -fpic to make the python shared libraries work, and I wanted to put
# the C++ gunk in subdirectories.  As a result, you'll need to
# tweak the CCFLAGS a bit, probably.

CC      = g++
CCFLAGS = -fpic -g -Wall -fmessage-length=0 -I/usr/local/include/python2.3
CCFLAGS = -fPIC -g -Wall -fmessage-length=0 -I/home/bwalenz/linux64/include/python2.4
CCFLAGS = -fpic -g -Wall -fmessage-length=0 -I/home/bwalenz/linux/include/python2.4
CCFLAGS = -fpic -g -Wall                    -I/usr/local/include/python2.3

SOURCE  = localalign/GF_ALN_dpaligner.cc \
          localalign/GF_ALN_local.cc \
          localalign/GF_ALN_overlap.cc \
          localalign/GF_ALN_loverlapper.cc \
          localalign/GF_ALN_pieceOlap.cc \
          localalign/localAlignerInterfacemodule.cc \
          halign/halign.cc \
          halign/halignmodule.cc

PY_EXES = python/AtacDriver.py \
          python/AtacFile.py \
          python/DNA.py \
          python/IdxStore.py \
          python/MatchRecord.py \
          python/MyFile.py \
          python/PerfectRuns.py \
          python/TrimMatchOverlaps.py \
          python/UniqueFilter.py \
          python/dedashMatches.py \
          python/fillIntraRunGaps.py \
          python/mkstats.py \
          python/squeezeIntraRunGaps.py

all:     localAlignerInterfacemodule.so halignmodule.so

install: localAlignerInterfacemodule.so halignmodule.so

clean:
	rm -f */*.o *.so python/*.pyc
real-clean:
	rm -f */*.o *.so python/*.pyc

localAlignerInterfacemodule.so:   localalign/GF_ALN_dpaligner.o \
                                  localalign/GF_ALN_local.o \
                                  localalign/GF_ALN_overlap.o \
                                  localalign/GF_ALN_loverlapper.o \
                                  localalign/GF_ALN_pieceOlap.o \
                                  localalign/localAlignerInterfacemodule.o
	$(CC) -shared -o $@ $^

halignmodule.so: halign/halignmodule.o halign/halign.o
	$(CC) -shared -o $@ $^

%.o: %.cc
	$(CC) $(CCFLAGS) -o $@ -c $<
