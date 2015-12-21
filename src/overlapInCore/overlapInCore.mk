#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := overlapInCore
SOURCES  := overlapInCore.C \
            overlapInCore-Build_Hash_Index.C \
            overlapInCore-Find_Overlaps.C \
            overlapInCore-Output.C \
            overlapInCore-Process_Overlaps.C \
            overlapInCore-Process_String_Overlaps.C

SRC_INCDIRS  := .. ../AS_UTL ../stores liboverlap

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
