
# The standard build system doesn't work well here.  We need to
# include -fpic to make the python shared libraries work, we need to
# include the python header, and I wanted to put the C++ gunk in
# subdirectories.  As a result, you'll need to tweak the CCFLAGS a
# bit, probably.

include ../../Make.compilers

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
	$(CC) -fpic $(CFLAGS_COMPILE) $(CFLAGS_PYTHON) -o $@ -c $<
