MODULE       :=    canu
TARGET       := libcanu.a
SOURCES      := utility/src/utility/runtime.C \
                \
                utility/src/utility/align-ssw.C \
                utility/src/utility/align-ssw-driver.C \
                utility/src/utility/edlib.C \
                \
                utility/src/utility/files.C \
                utility/src/utility/files-buffered.C \
                utility/src/utility/files-compressed.C \
                utility/src/utility/files-memoryMapped.C \
                \
                utility/src/utility/logging.C \
                \
                utility/src/utility/strings.C \
                \
                utility/src/utility/system.C \
                utility/src/utility/system-stackTrace.C \
                \
                utility/src/utility/sequence.C \
                \
                utility/src/utility/types.C \
                \
                utility/src/utility/kmers-exact.C \
                utility/src/utility/kmers-files.C \
                utility/src/utility/kmers-histogram.C \
                utility/src/utility/kmers-reader.C \
                utility/src/utility/kmers-writer-block.C \
                utility/src/utility/kmers-writer-stream.C \
                utility/src/utility/kmers-writer.C \
                utility/src/utility/kmers.C \
                \
                utility/src/utility/bits.C \
                \
                utility/src/utility/hexDump.C \
                utility/src/utility/md5.C \
                utility/src/utility/mt19937ar.C \
                utility/src/utility/speedCounter.C \
                utility/src/utility/sweatShop.C \
                \
                correction/computeGlobalScore.C \
                correction/falconConsensus.C \
                correction/falconConsensus-alignTag.C \
                \
                stores/sqCache.C \
                stores/sqLibrary.C \
                stores/sqReadData.C \
                stores/sqReadDataWriter.C \
                stores/sqStore.C \
                stores/sqStoreBlob.C \
                stores/sqStoreConstructor.C \
                stores/sqStoreInfo.C \
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
                stores/objectStore.C \
                \
                overlapInCore/overlapReadCache.C \
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
                utgcns/libpbutgcns/AlnGraphBoost.C  \
                \
                gfa/gfa.C \
                gfa/bed.C


ifeq (${BUILDSTACKTRACE}, 1)
SOURCES      += utility/src/utility/libbacktrace/atomic.c \
                utility/src/utility/libbacktrace/backtrace.c \
                utility/src/utility/libbacktrace/dwarf.c \
                utility/src/utility/libbacktrace/elf.c \
                utility/src/utility/libbacktrace/fileline.c \
                utility/src/utility/libbacktrace/mmap.c \
                utility/src/utility/libbacktrace/mmapio.c \
                utility/src/utility/libbacktrace/posix.c \
                utility/src/utility/libbacktrace/print.c \
                utility/src/utility/libbacktrace/simple.c \
                utility/src/utility/libbacktrace/sort.c \
                utility/src/utility/libbacktrace/state.c \
                utility/src/utility/libbacktrace/unknown.c
endif



SRC_INCDIRS  := . \
                utility/src/utility \
                stores \
                stores/libsnappy \
                alignment \
                utgcns/libNDalign \
                utgcns/libcns \
                utgcns/libpbutgcns \
                utgcns/libNDFalcon \
                utgcns/libboost \
                overlapInCore \
                overlapInCore/liboverlap

SUBMAKEFILES := stores/dumpBlob.mk \
                stores/ovStoreBuild.mk \
                stores/ovStoreConfig.mk \
                stores/ovStoreBucketizer.mk \
                stores/ovStoreSorter.mk \
                stores/ovStoreIndexer.mk \
                stores/ovStoreDump.mk \
                stores/ovStoreStats.mk \
                stores/sqStoreCreate.mk \
                stores/sqStoreDumpFASTQ.mk \
                stores/sqStoreDumpMetaData.mk \
                stores/tgStoreCompress.mk \
                stores/tgStoreDump.mk \
                stores/tgStoreLoad.mk \
                stores/tgTigDisplay.mk \
                stores/loadCorrectedReads.mk \
                stores/loadTrimmedReads.mk \
                stores/loadErates.mk \
                \
                meryl/src/meryl/meryl.mk \
                \
                overlapInCore/overlapInCore.mk \
                overlapInCore/overlapInCorePartition.mk \
                overlapInCore/overlapConvert.mk \
                overlapInCore/overlapImport.mk \
                overlapInCore/overlapPair.mk \
                overlapInCore/edalign.mk \
                \
                overlapInCore/liboverlap/prefixEditDistance-matchLimitGenerate.mk \
                \
                mhap/mhapConvert.mk \
                \
                minimap/mmapConvert.mk \
                \
                correction/filterCorrectionOverlaps.mk \
                correction/generateCorrectionLayouts.mk \
                correction/filterCorrectionLayouts.mk \
                correction/falconsense.mk \
                \
                correction/errorEstimate.mk \
                \
                haplotyping/splitHaplotype.mk \
                \
                overlapBasedTrimming/trimReads.mk \
                overlapBasedTrimming/splitReads.mk \
                overlapBasedTrimming/mergeRanges.mk \
                \
                overlapAlign/overlapAlign.mk \
                \
                overlapErrorAdjustment/findErrors.mk \
                overlapErrorAdjustment/fixErrors.mk \
                overlapErrorAdjustment/findErrors-Dump.mk \
                overlapErrorAdjustment/correctOverlaps.mk \
                \
                bogart/bogart.mk \
                bogart/layoutReads.mk \
                \
                bogus/bogus.mk \
                \
                utgcns/utgcns.mk \
                \
                gfa/alignGFA.mk
