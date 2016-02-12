#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := overlapPair
SOURCES  := overlapPair.C

SRC_INCDIRS  := .. ../AS_UTL ../stores ../meryl/libleaff liboverlap ../utgcns/libNDFalcon

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lleaff -lcanu
TGT_PREREQS := libleaff.a libcanu.a

SUBMAKEFILES :=
