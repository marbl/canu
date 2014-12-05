
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

TARGET   := libkmer.a

SOURCES  := existDB-create-from-fasta.C \
            existDB-create-from-meryl.C \
            existDB-create-from-sequence.C \
            existDB-state.C \
            existDB.C \
            positionDB-access.C \
            positionDB-dump.C \
            positionDB-file.C \
            positionDB-mismatch.C \
            positionDB-sort.C \
            positionDB.C

SRC_INCDIRS  := ../libutil ../libbio ../libseq ../libmeryl

SUBMAKEFILES := existDB.mk posDB.mk kmer-mask.mk percentCovered.mk
