
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

TARGET   := meryl
SOURCES  := args.C \
            binaryOp.C \
            build-threads.C \
            build.C \
            dump.C \
            estimate.C \
            merge.C \
            meryl.C \
            unaryOp.C

#           compare-counts.C
#           mervin.C
#           simple.C

SRC_INCDIRS := ../libutil ../libbio ../libseq ../libkmer ../libmeryl

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lseq -lbio -lkmer -lmeryl -lutil
TGT_PREREQS := libutil.a libbio.a libkmer.a libmeryl.a

#           merge.qsort.C \
#           merge.listmerge.C \

#           asmMerQC-regions.C \
#           asmMerQC.C \
#           testPositionBias.C \
#           m.C \
#           mapMers-depth.C \
#           mapMers.C \
#           maskMers.C \

SUBMAKEFILES :=

