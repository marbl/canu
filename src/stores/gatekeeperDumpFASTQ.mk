
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := gatekeeperDumpFASTQ
SOURCES  := gatekeeperDumpFASTQ.C

SRC_INCDIRS := .. ../stores ../AS_UTL ../AS_OBT

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lCA
TGT_PREREQS := libCA.a

SUBMAKEFILES :=
