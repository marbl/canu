
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
            alignment/overlap.C \
            \
            overlapErrorAdjustment/analyzeAlignment.C \
            \
            overlapInCore/liboverlap/Binomial_Bound.C \
            overlapInCore/liboverlap/prefixEditDistance.C \
            overlapInCore/liboverlap/prefixEditDistance-allocateMoreSpace.C \
            overlapInCore/liboverlap/prefixEditDistance-extend.C \
            overlapInCore/liboverlap/prefixEditDistance-forward.C \
            overlapInCore/liboverlap/prefixEditDistance-reverse.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-5000.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4900.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4800.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4700.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4600.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4500.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4400.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4300.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4200.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4100.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-4000.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3900.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3800.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3700.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3600.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3500.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3400.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3300.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3200.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3100.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-3000.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2900.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2800.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2700.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2600.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2500.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2400.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2300.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2200.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2100.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-2000.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1900.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1800.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1700.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1600.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1500.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1400.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1300.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1200.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1100.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-1000.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0900.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0800.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0700.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0600.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0500.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0400.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0300.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0200.C \
            overlapInCore/liboverlap/prefixEditDistance-matchLimit-0100.C \
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

#  The last two are for aligners/overlap.C
SRC_INCDIRS  := . AS_UTL stores alignment overlapInCore/liboverlap utgcns/libcns   meryl/libleaff overlapInCore

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
                meryl/libleaff.mk \
                meryl/leaff.mk \
                meryl/meryl.mk \
                meryl/simple.mk \
                meryl/estimate-mer-threshold.mk \
                \
                overlapInCore/overlapInCore.mk \
                overlapInCore/overlapInCorePartition.mk \
                overlapInCore/overlapConvert.mk \
                overlapInCore/overlapImport.mk \
                overlapInCore/overlapPair.mk \
                \
                overlapInCore/liboverlap/prefixEditDistance-matchLimitGenerate.mk \
                \
                mhap/mhap.mk \
                mhap/mhapConvert.mk \
                \
                correction/filterCorrectionOverlaps.mk \
                correction/generateCorrectionLayouts.mk \
                \
                falcon_sense/createFalconSenseInputs.mk \
                \
                overlapBasedTrimming/trimReads.mk \
                overlapBasedTrimming/splitReads.mk \
                \
                overlapErrorAdjustment/findErrors.mk \
                overlapErrorAdjustment/correctOverlaps.mk \
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
