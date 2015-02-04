
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := libleaff.a
SOURCES  := libleaff/fastaFile.C \
            libleaff/fastaStdin.C \
            libleaff/fastqFile.C \
            libleaff/fastqStdin.C \
            libleaff/gkStoreFile.C \
            libleaff/merStream.C \
            libleaff/seqCache.C \
            libleaff/seqFactory.C \
            libleaff/seqStore.C \
            libleaff/seqStream.C \
            libleaff/sffFile.C

#            libleaff/selftest.C
#            libleaff/test-merStream.C
#            libleaff/test-seqCache.C
#            libleaff/test-seqStream.C

SRC_INCDIRS  := .. ../AS_UTL ../stores

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := 
TGT_PREREQS := 

SUBMAKEFILES :=

