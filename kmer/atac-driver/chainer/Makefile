
# Ross's build system wasn't doing well here.  We need to include
# -fpic to make the python shared libraries work, and I wanted to put
# the C++ gunk in subdirectories.

CC      = g++
CCFLAGS = -fpic -O -Wall -I/usr/local/include/python2.3

SOURCE  = localalign/GF_ALN_dpaligner.cc \
          localalign/GF_ALN_local.cc \
          localalign/GF_ALN_overlap.cc \
          localalign/GF_ALN_qvaligner.cc \
          localalign/GF_ALN_loverlapper.cc \
          localalign/GF_ALN_pieceOlap.cc \
          localalign/localAlignerInterface.cc \
          localalign/localAlignerInterfacemodule.cc \
          halign/halign.cc \
          halign/halignmodule.cc

all: localAlignerInterfacemodule.so halignmodule.so

clean:
	rm -f */*.o *.so *.pyc
real-clean:
	rm -f */*.o *.so *.pyc

localAlignerInterfacemodule.so:   localalign/GF_ALN_dpaligner.o \
                                  localalign/GF_ALN_local.o \
                                  localalign/GF_ALN_overlap.o \
                                  localalign/GF_ALN_qvaligner.o \
                                  localalign/GF_ALN_loverlapper.o \
                                  localalign/GF_ALN_pieceOlap.o \
                                  localalign/localAlignerInterface.o \
                                  localalign/localAlignerInterfacemodule.o
	$(CC) -shared -o $@ $^

halignmodule.so: halign/halignmodule.o halign/halign.o
	$(CC) -shared -o $@ $^

%.o: %.cc
	${CC} ${CCFLAGS} -o $@ -c $<

%.so:
	ar 