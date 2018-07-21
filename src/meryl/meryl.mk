
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)
endif

TARGET   := meryl
SOURCES  := meryl.C \
            merylCountArray.C \
            merylInput.C \
            merylOp-count.C \
            merylOp-countSimple.C \
            merylOp-histogram.C \
            merylOp-nextMer.C \
            merylOp.C

SRC_INCDIRS  := . .. ../utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
