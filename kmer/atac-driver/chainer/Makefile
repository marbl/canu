
# Ross's build system wasn't doing well here.  We need to include
# -fpic to make the python shared libraries work, and I wanted to put
# the C++ gunk in subdirectories.  As a result, you'll need to
# tweak the CCFLAGS a bit, probably.

CC      = g++
CCFLAGS = -fpic -g -Wall -fmessage-length=0 -I/usr/local/include/python2.3
CCFLAGS = -fPIC -g -Wall -fmessage-length=0 -I/home/bwalenz/linux64/include/python2.4
CCFLAGS = -fpic -g -Wall -fmessage-length=0 -I/home/bwalenz/linux/include/python2.4

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

all: localAlignerInterfacemodule.so halignmodule.so atacdriver.sh

clean:
	rm -f */*.o *.so python/*.pyc atacdriver.sh
real-clean:
	rm -f */*.o *.so python/*.pyc atacdriver.sh

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

atacdriver.sh: Makefile
	echo \#!/bin/sh > atacdriver.sh
	echo export PYTHONPATH=`pwd`:`pwd`/python >> atacdriver.sh
	echo python `pwd`/python/AtacDriver.py $$\* >> atacdriver.sh


%.o: %.cc
	$(CC) $(CCFLAGS) -o $@ -c $<

%.so:
	ar 
