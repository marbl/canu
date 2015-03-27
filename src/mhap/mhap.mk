#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mhap-0.1-ob.jar
SOURCES  := mhap-0.1-ob.jar lib/guava-16.0.jar lib/jaligner.jar

SRC_INCDIRS  := 

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lCA
TGT_PREREQS := libCA.a

SUBMAKEFILES :=
