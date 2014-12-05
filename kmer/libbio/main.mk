
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

#CXX =g++46
#CXXFLAGS := -fopenmp -pthread -D_GLIBCXX_PARALLEL -O3 -fPIC -m64 -pipe -Wno-write-strings -Wall -Wno-char-subscripts -mtune=native -march=native -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
#LDFLAGS  := -fopenmp -pthread -lm -lthr

TARGET   := libbio.a

SOURCES  := alphabet.c \
            alphabet-acgtspace.c \
            alphabet-colorspace.c \
            halign.c \
            kmer.C \
            reversecomplement.c

SRC_INCDIRS = ../libutil

SUBMAKEFILES := 
