
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := meryl
SOURCES  := meryl-args.C \
            meryl-binaryOp.C \
            meryl-build.C \
            meryl-build-threads.C \
            meryl-dump.C \
            meryl-estimate.C \
            meryl-merge.C \
            meryl-unaryOp.C \
            meryl.C

SRC_INCDIRS  := .. ../AS_UTL libleaff

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lleaff -lcanu
TGT_PREREQS := libleaff.a libcanu.a

SUBMAKEFILES :=

