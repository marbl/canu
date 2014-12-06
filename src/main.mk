
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := libCA.a

SOURCES  := AS_global.C \
            \
            AS_UTL/AS_UTL_GPL.C \
            AS_UTL/AS_UTL_Hash.C \
            AS_UTL/AS_UTL_Var.C \
            AS_UTL/AS_UTL_alloc.C \
            AS_UTL/AS_UTL_decodeRange.C \
            AS_UTL/AS_UTL_fasta.C \
            AS_UTL/AS_UTL_fileIO.C \
            AS_UTL/AS_UTL_heap.C \
            AS_UTL/AS_UTL_histo.C \
            AS_UTL/AS_UTL_interval.C \
            AS_UTL/AS_UTL_rand.C \
            AS_UTL/AS_UTL_reverseComplement.C \
            AS_UTL/AS_UTL_stackTrace.C \
            AS_UTL/AS_UTL_testhisto.C \
            \
            stores/gkStore.C \
            \
            meryl/libmeryl.C

#  Not complete, just snippets of code
#            AS_UTL/decodeBooleanString.C
#            AS_UTL/AS_UTL_stackTraceTest.C
#
#  Test code
#            AS_UTL/testHashTable.C
#            AS_UTL/testRand.C
#            AS_UTL/testVar.C
#            AS_UTL/memoryMappedFileTest.C

SRC_INCDIRS  := . AS_UTL

#  Cannot be the 'main.mk' from subdirectories, as that redefines libCA.a to
#  update just those objects.  Instead, this must list each binary indifidually.
SUBMAKEFILES := stores/gatekeeperCreate.mk \
                stores/gatekeeperDumpFASTQ.mk \
                meryl/main.mk
