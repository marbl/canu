
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := leaff
SOURCES  := leaff-blocks.C \
            leaff-duplicates.C \
            leaff-gc.C \
            leaff-partition.C \
            leaff-simulate.C \
            leaff-statistics.C \
            leaff.C

SRC_INCDIRS  := .. ../AS_UTL libleaff

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lleaff -lcanu
TGT_PREREQS := libleaff.a libcanu.a

SUBMAKEFILES :=

