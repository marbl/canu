
# The standard build system doesn't work well here.  We need to
# include -fpic to make the python shared libraries work, we need to
# include the python header, and I wanted to put the C++ gunk in
# subdirectories.  As a result, you might need to tweak the CCFLAGS a
# bit.

include ../../Make.compilers

SOURCE  = localalign/GF_ALN_dpaligner.C \
          localalign/GF_ALN_local.C \
          localalign/GF_ALN_overlap.C \
          localalign/GF_ALN_loverlapper.C \
          localalign/GF_ALN_pieceOlap.C \
          localalign/localAlignerInterfacemodule.C \
          halign/halign.C \
          halign/halignmodule.C

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
	$(CXX) $(CXXSHARED) -o $@ $^

halignmodule.so: halign/halignmodule.o halign/halign.o
	$(CXX) $(CXXSHARED) -o $@ $^

%.o: %.C
	$(CXX) -fPIC $(CFLAGS_COMPILE) $(CFLAGS_PYTHON) -o $@ -c $<
