
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.

ifeq "$(strip ${DESTDIR})" ""
  DESTDIR      := 
endif

ifeq "$(strip ${PREFIX})" ""
  ifeq "$(strip ${DESTDIR})" ""
    PREFIX     := $(realpath ..)
  else
    PREFIX     := /canu
  endif
endif

ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := $(DESTDIR)$(PREFIX)/$(OSTYPE)-$(MACHINETYPE)/obj
endif

ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := $(DESTDIR)$(PREFIX)/$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libcanu.a

SOURCES      := AS_global.C \
                \
                AS_UTL/AS_UTL_decodeRange.C \
                AS_UTL/AS_UTL_fasta.C \
                AS_UTL/AS_UTL_fileIO.C \
                AS_UTL/AS_UTL_reverseComplement.C \
                AS_UTL/AS_UTL_stackTrace.C \
                \
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
                stores/gkStoreEncode.C \
                \
                stores/ovOverlap.C \
                stores/ovStore.C \
                stores/ovStoreFile.C \
                \
                stores/tgStore.C \
                stores/tgTig.C \
                stores/tgTigSizeAnalysis.C \
                stores/tgTigMultiAlignDisplay.C \
                \
                meryl/libmeryl.C \
                \
                overlapInCore/overlapReadCache.C \
                \
                overlapErrorAdjustment/analyzeAlignment.C \
                \
                overlapInCore/liboverlap/Binomial_Bound.C \
                overlapInCore/liboverlap/Display_Alignment.C \
                overlapInCore/liboverlap/prefixEditDistance.C \
                overlapInCore/liboverlap/prefixEditDistance-allocateMoreSpace.C \
                overlapInCore/liboverlap/prefixEditDistance-extend.C \
                overlapInCore/liboverlap/prefixEditDistance-forward.C \
                overlapInCore/liboverlap/prefixEditDistance-reverse.C \
                \
                utgcns/libNDalign/NDalign.C \
                \
                utgcns/libNDalign/Binomial_Bound.C \
                utgcns/libNDalign/NDalgorithm.C \
                utgcns/libNDalign/NDalgorithm-allocateMoreSpace.C \
                utgcns/libNDalign/NDalgorithm-extend.C \
                utgcns/libNDalign/NDalgorithm-forward.C \
                utgcns/libNDalign/NDalgorithm-reverse.C \
                \
                utgcns/libcns/abAbacus-addRead.C \
                utgcns/libcns/abAbacus-appendBases.C \
                utgcns/libcns/abAbacus-applyAlignment.C \
                utgcns/libcns/abAbacus-baseCall.C \
                utgcns/libcns/abAbacus-mergeRefine.C \
                utgcns/libcns/abAbacus-refine.C \
                utgcns/libcns/abAbacus-refreshMultiAlign.C \
                utgcns/libcns/abAbacus.C \
                utgcns/libcns/abColumn.C \
                utgcns/libcns/abMultiAlign.C \
                utgcns/libcns/unitigConsensus.C \
                utgcns/libpbutgcns/Alignment.C	\
                utgcns/libpbutgcns/AlnGraphBoost.C  \
                utgcns/libpbutgcns/SimpleAligner.C \
                utgcns/libNDFalcon/dw.C \
                \
                falcon_sense/libfalcon/falcon.C \
                falcon_sense/libfalcon/kmer_lookup.C

SRC_INCDIRS  := . \
                AS_UTL \
                stores \
                alignment \
                utgcns/libNDalign \
                utgcns/libcns \
                utgcns/libpbutgcns \
                utgcns/libNDFalcon \
                utgcns/libboost \
                meryl/libleaff \
                overlapInCore \
                overlapInCore/liboverlap \
                falcon_sense/libfalcon

SUBMAKEFILES := stores/gatekeeperCreate.mk \
                stores/gatekeeperDumpFASTQ.mk \
                stores/gatekeeperDumpMetaData.mk \
                stores/gatekeeperPartition.mk \
                stores/ovStoreBuild.mk \
                stores/ovStoreBucketizer.mk \
                stores/ovStoreSorter.mk \
                stores/ovStoreIndexer.mk \
                stores/ovStoreDump.mk \
                stores/ovStoreStats.mk \
                stores/tgStoreCompress.mk \
                stores/tgStoreDump.mk \
                stores/tgStoreLoad.mk \
                stores/tgStoreFilter.mk \
                stores/tgStoreCoverageStat.mk \
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
                minimap/mmapConvert.mk \
                \
                correction/filterCorrectionOverlaps.mk \
                correction/generateCorrectionLayouts.mk \
                correction/readConsensus.mk \
                correction/errorEstimate.mk \
                \
                falcon_sense/createFalconSenseInputs.mk \
                falcon_sense/falcon_sense.mk \
                \
                overlapBasedTrimming/trimReads.mk \
                overlapBasedTrimming/splitReads.mk \
                \
                overlapErrorAdjustment/findErrors.mk \
                overlapErrorAdjustment/findErrors-Dump.mk \
                overlapErrorAdjustment/correctOverlaps.mk \
                \
                bogart/bogart.mk \
                bogart/buildGraph.mk \
                \
                bogus/bogus.mk \
                \
                erateEstimate/erateEstimate.mk \
                \
                utgcns/utgcns.mk \
                \
                fastq-utilities/fastqAnalyze.mk \
                fastq-utilities/fastqSample.mk \
                fastq-utilities/fastqSimulate.mk \
                fastq-utilities/fastqSimulate-sort.mk
