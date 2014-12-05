
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

#CXX = g++46
#CC  = gcc46
#CXXFLAGS := -fopenmp -D_GLIBCXX_PARALLEL -O3 -fPIC -m64 -pipe -Wno-write-strings
#LDFLAGS  := -fopenmp -lm

TARGET   := libutil.a

#bigQueue.C

SOURCES  := bitPackedArray.C \
            bitPackedFile.C \
            fibonacciNumbers.C \
            readBuffer.C \
            recordFile.C \
            speedCounter.C \
            sweatShop.C \
            file.c \
            md5.c \
            palloc.c \
            qsort_mt.c \
            util.c \
            mt19937ar/mt19937ar.c

#           unaryEncodingTester.C
#           bzipBuffer.C

SRC_INCDIRS  :=
SUBMAKEFILES := 
