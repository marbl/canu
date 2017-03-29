
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
                AS_UTL/hexDump.C \
                AS_UTL/md5.C \
                AS_UTL/mt19937ar.C \
                AS_UTL/readBuffer.C \
                AS_UTL/speedCounter.C \
                AS_UTL/sweatShop.C \
                AS_UTL/timeAndSize.C \
                AS_UTL/kMer.C \
                \
                falcon_sense/libfalcon/falcon.C \
                \
                stores/gkStore.C \
                stores/gkStoreEncode.C \
                \
                stores/ovOverlap.C \
                stores/ovStore.C \
                stores/ovStoreWriter.C \
                stores/ovStoreFilter.C \
                stores/ovStoreFile.C \
                stores/ovStoreHistogram.C \
                \
                stores/tgStore.C \
                stores/tgTig.C \
                stores/tgTigSizeAnalysis.C \
                stores/tgTigMultiAlignDisplay.C \
                \
                stores/libsnappy/snappy-sinksource.cc \
                stores/libsnappy/snappy-stubs-internal.cc \
                stores/libsnappy/snappy.cc \
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
                overlapInCore/libedlib/edlib.C \
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
                utgcns/libNDFalcon/dw.C \
                \
                meryl/libkmer/existDB-create-from-fasta.C \
                meryl/libkmer/existDB-create-from-meryl.C \
                meryl/libkmer/existDB-create-from-sequence.C \
                meryl/libkmer/existDB-state.C \
                meryl/libkmer/existDB.C \
                meryl/libkmer/positionDB-access.C \
                meryl/libkmer/positionDB-dump.C \
                meryl/libkmer/positionDB-file.C \
                meryl/libkmer/positionDB-mismatch.C \
                meryl/libkmer/positionDB-sort.C \
                meryl/libkmer/positionDB.C



ifeq (${BUILDSTACKTRACE}, 1)
SOURCES      += AS_UTL/libbacktrace/atomic.c \
                AS_UTL/libbacktrace/backtrace.c \
                AS_UTL/libbacktrace/dwarf.c \
                AS_UTL/libbacktrace/elf.c \
                AS_UTL/libbacktrace/fileline.c \
                AS_UTL/libbacktrace/mmap.c \
                AS_UTL/libbacktrace/mmapio.c \
                AS_UTL/libbacktrace/posix.c \
                AS_UTL/libbacktrace/print.c \
                AS_UTL/libbacktrace/simple.c \
                AS_UTL/libbacktrace/sort.c \
                AS_UTL/libbacktrace/state.c \
                AS_UTL/libbacktrace/unknown.c
endif



SRC_INCDIRS  := . \
                AS_UTL \
                stores \
                stores/libsnappy \
                alignment \
                utgcns/libNDalign \
                utgcns/libcns \
                utgcns/libpbutgcns \
                utgcns/libNDFalcon \
                utgcns/libboost \
                meryl/libleaff \
                overlapInCore \
                overlapInCore/libedlib \
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
                meryl/maskMers.mk \
                meryl/simple.mk \
                meryl/estimate-mer-threshold.mk \
                meryl/existDB.mk \
                meryl/positionDB.mk \
                \
                merTrim/merTrim.mk \
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
                \
                bogus/bogus.mk \
                \
                erateEstimate/erateEstimate.mk \
                \
                utgcns/utgcns.mk \
                \
                utgcns/alignGFA.mk \
                \
                fastq-utilities/fastqAnalyze.mk \
                fastq-utilities/fastqSample.mk \
                fastq-utilities/fastqSimulate.mk \
                fastq-utilities/fastqSimulate-sort.mk
