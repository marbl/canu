
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := libCA.a

SOURCES  := AS_UTL_GPL.C \
            AS_UTL_Hash.C \
            AS_UTL_Var.C \
            AS_UTL_alloc.C \
            AS_UTL_decodeRange.C \
            AS_UTL_fasta.C \
            AS_UTL_fileIO.C \
            AS_UTL_heap.C \
            AS_UTL_histo.C \
            AS_UTL_interval.C \
            AS_UTL_rand.C \
            AS_UTL_reverseComplement.C \
            AS_UTL_stackTrace.C \
            kMer.C \

#            decodeBooleanString.C \
#            AS_UTL_stackTraceTest.C \
#            AS_UTL_testhisto.C \
#            memoryMappedFileTest.C \
#            testHashTable.C \
#            testRand.C \
#            testVar.C \

SRC_INCDIRS  := . ..
SUBMAKEFILES := 
