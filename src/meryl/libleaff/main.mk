
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

#CXX = g++46
#CXXFLAGS := -fopenmp -D_GLIBCXX_PARALLEL -O3 -fPIC -m64 -pipe -Wno-write-strings
#LDFLAGS  := -fopenmp -lm

TARGET   := libleaff.a
SOURCES  := fastaFile.C \
            fastaStdin.C \
            fastqFile.C \
            fastqStdin.C \
            gkStoreFile.C \
            merStream.C \
            seqCache.C \
            seqFactory.C \
            seqStore.C \
            seqStream.C \
            sffFile.C

SRC_INCDIRS  := ../.. ../../AS_UTL ../../stores

#SUBMAKEFILES := test-merStream.mk test-seqCache.mk test-seqStream.mk


