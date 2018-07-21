
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)
endif

TARGET   := libleaff.a
SOURCES  := libleaff/dnaAlphabets.C \
            libleaff/fastaFile.C \
            libleaff/fastaStdin.C \
            libleaff/fastqFile.C \
            libleaff/fastqStdin.C \
            libleaff/sqStoreFile.C \
            libleaff/merStream.C \
            libleaff/seqCache.C \
            libleaff/seqFactory.C \
            libleaff/seqStore.C \
            libleaff/seqStream.C \
            libleaff/sffFile.C \
            utility/bitEncodings.C \
            utility/bitPackedFile.C \
            utility/bitPackedArray.C \
            utility/kMer.C \
            libkmer/existDB-create-from-fasta.C \
            libkmer/existDB-create-from-meryl.C \
            libkmer/existDB-create-from-sequence.C \
            libkmer/existDB-state.C \
            libkmer/existDB.C \
            libkmer/positionDB-access.C \
            libkmer/positionDB-dump.C \
            libkmer/positionDB-file.C \
            libkmer/positionDB-mismatch.C \
            libkmer/positionDB-sort.C \
            libkmer/positionDB.C \
            libmeryl-san.C

#            libleaff/selftest.C
#            libleaff/test-merStream.C
#            libleaff/test-seqCache.C
#            libleaff/test-seqStream.C

SRC_INCDIRS  := . ./libleaff ./utility ./libkmer .. ../utility ../stores

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := 
TGT_PREREQS := 

SUBMAKEFILES :=

