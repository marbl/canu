
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := libCA.a

#  unitigger uses AS_UTL_alloc, AS_UTL_Var and AS_UTL_Hash
#  AS_UTL_Hash needs AS_UTL_heap

SOURCES  := AS_global.C \
            \
            AS_UTL/AS_UTL_decodeRange.C \
            AS_UTL/AS_UTL_fasta.C \
            AS_UTL/AS_UTL_fileIO.C \
            AS_UTL/AS_UTL_reverseComplement.C \
            AS_UTL/AS_UTL_stackTrace.C \
            \
            AS_UTL/AS_UTL_Var.C \
            AS_UTL/AS_UTL_Hash.C \
            AS_UTL/AS_UTL_heap.C \
            AS_UTL/AS_UTL_alloc.C \
            \
            AS_UTL/bitEncodings.C \
            AS_UTL/bitPackedFile.C \
            AS_UTL/bitPackedArray.C \
            AS_UTL/dnaAlphabets.C \
            AS_UTL/md5.C \
            AS_UTL/mt19937ar.C \
            AS_UTL/readBuffer.C \
            AS_UTL/speedCounter.C \
            AS_UTL/sweatShop.C \
            AS_UTL/timeAndSize.C \
            AS_UTL/kMer.C \
            \
            stores/gkStore.C \
            \
            stores/ovOverlap.C \
            stores/ovStore.C \
            stores/ovStoreFile.C \
            \
            stores/tgStore.C \
            stores/tgTig.C \
            stores/tgTigSizeAnalysis.C \
            stores/tgTigMultiAlignment.C \
            stores/tgTigMultiAlignDisplay.C \
            \
            meryl/libmeryl.C \
            \
            alignment/alignment-drivers.C \
            alignment/analyze-affine-overlap.C \
            alignment/brute-force-dp.C \
            alignment/dp-align.C \
            alignment/find-local-overlap.C \
            alignment/find-local-segments.C \
            alignment/fix-overlapping-pieces.C \
            alignment/local-overlapper.C \
            \
            overlapInCore/liboverlap/Binomial_Bound.C \
            \
            utgcns/libcns/abAbacus-addRead.C \
            utgcns/libcns/abAbacus-applyAlignment.C \
            utgcns/libcns/abAbacus-baseCall.C \
            utgcns/libcns/abAbacus-mergeRefine.C \
            utgcns/libcns/abAbacus-refine.C \
            utgcns/libcns/abAbacus-refreshMultiAlign.C \
            utgcns/libcns/abAbacus.C \
            utgcns/libcns/abColumn.C \
            utgcns/libcns/abMultiAlign.C \
            utgcns/libcns/unitigConsensus.C


#  Not complete, just snippets of code
#            AS_UTL/decodeBooleanString.C
#            AS_UTL/AS_UTL_stackTraceTest.C
#
#  Test code
#            AS_UTL/testHashTable.C
#            AS_UTL/testRand.C
#            AS_UTL/testVar.C
#            AS_UTL/memoryMappedFileTest.C

SRC_INCDIRS  := . AS_UTL stores alignment overlapInCore/liboverlap utgcns/libcns

#  Cannot be the 'main.mk' from subdirectories, as that redefines libCA.a to
#  update just those objects.  Instead, this must list each binary indifidually.
SUBMAKEFILES := stores/gatekeeperCreate.mk \
                stores/gatekeeperDumpFASTQ.mk \
                stores/gatekeeperDumpMetaData.mk \
                stores/gatekeeperPartition.mk \
                stores/ovStoreBuild.mk \
                stores/ovStoreBucketizer.mk \
                stores/ovStoreSorter.mk \
                stores/ovStoreIndexer.mk \
                stores/ovStoreDump.mk \
                stores/tgStoreDump.mk \
                stores/tgTigDisplay.mk \
                \
                meryl/main.mk \
                \
                overlapInCore/overlapInCore.mk \
                overlapInCore/overlapInCorePartition.mk \
                overlapInCore/overlapConvert.mk \
                overlapInCore/overlapPair.mk \
                \
                AS_OBT/deduplicate.mk \
                AS_OBT/finalTrim.mk \
                AS_OBT/chimera.mk \
                \
                AS_OVL/findErrors.mk \
                AS_OVL/correctOverlaps.mk \
                \
                bogart/bogart.mk \
                \
                bogus/bogus.mk \
                \
                erateEstimate/erateEstimate.mk \
                \
                unitigger/unitigger.mk \
                \
                utgcns/utgcns.mk \
                \
                fastq-utilities/fastqAnalyze.mk \
                fastq-utilities/fastqSample.mk \
                fastq-utilities/fastqSimulate.mk \
                fastq-utilities/fastqSimulate-sort.mk
