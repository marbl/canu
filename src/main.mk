MODULE       :=    canu
TARGET       := libcanu.a
SOURCES      := utility/src/align/align-ksw2-driver.C \
                utility/src/align/align-ksw2-extz.C \
                utility/src/align/align-ksw2-extz2-sse.C \
                utility/src/align/align-parasail-driver.C \
                utility/src/align/align-ssw-driver.C \
                utility/src/align/align-ssw.C \
                utility/src/align/edlib.C \
                \
                utility/src/bits/fibonacci-v1.C \
                utility/src/bits/hexDump-v1.C \
                utility/src/bits/stuffedBits-v1-binary.C \
                utility/src/bits/stuffedBits-v1-bits.C \
                utility/src/bits/stuffedBits-v1-delta.C \
                utility/src/bits/stuffedBits-v1-gamma.C \
                utility/src/bits/stuffedBits-v1-golomb.C \
                utility/src/bits/stuffedBits-v1-omega.C \
                utility/src/bits/stuffedBits-v1-unary.C \
                utility/src/bits/stuffedBits-v1-zeckendorf.C \
                utility/src/bits/stuffedBits-v1.C \
                utility/src/bits/wordArray-v1.C \
                \
                utility/src/datastructures/keyAndValue-v1.C \
                utility/src/datastructures/splitToWords-v1.C \
                utility/src/datastructures/stringList-v1.C \
                utility/src/datastructures/strings-v1.C \
                utility/src/datastructures/types-v1.C \
                \
                utility/src/files/accessing-v1.C \
                utility/src/files/buffered-v1-reading.C \
                utility/src/files/buffered-v1-writing.C \
                utility/src/files/compressed-v1-reading.C \
                utility/src/files/compressed-v1-writing.C \
                utility/src/files/compressed-v1.C \
                utility/src/files/fasta-fastq-v1.C \
                utility/src/files/files-v1.C \
                utility/src/files/memoryMapped-v1.C \
                utility/src/files/readLine-v0.C \
                utility/src/files/readLine-v1.C \
                utility/src/files/reading-v1.C \
                utility/src/files/writing-v1.C \
                \
                utility/src/kmers-v1/kmers-exact.C \
                utility/src/kmers-v1/kmers-files.C \
                utility/src/kmers-v1/kmers-histogram.C \
                utility/src/kmers-v1/kmers-reader.C \
                utility/src/kmers-v1/kmers-writer-block.C \
                utility/src/kmers-v1/kmers-writer-stream.C \
                utility/src/kmers-v1/kmers-writer.C \
                utility/src/kmers-v1/kmers.C \
                \
                utility/src/kmers-v2/kmers-exact.C \
                utility/src/kmers-v2/kmers-files.C \
                utility/src/kmers-v2/kmers-histogram.C \
                utility/src/kmers-v2/kmers-reader-dump.C \
                utility/src/kmers-v2/kmers-reader.C \
                utility/src/kmers-v2/kmers-writer-block.C \
                utility/src/kmers-v2/kmers-writer-stream.C \
                utility/src/kmers-v2/kmers-writer.C \
                utility/src/kmers-v2/kmers.C \
                \
                utility/src/math/md5-v1.C \
                utility/src/math/mt19937ar-v1.C \
                utility/src/math/sampledDistribution-v1.C \
                \
                utility/src/parasail/cpuid.c \
                utility/src/parasail/memory.c \
                utility/src/parasail/sg.c \
                utility/src/parasail/sg_trace.c \
                utility/src/parasail/sg_qx_dispatch.c \
                utility/src/parasail/sg_qb_de_dispatch.c \
                utility/src/parasail/sg_qe_db_dispatch.c \
                utility/src/parasail/cigar.c \
                \
                utility/src/sequence/dnaSeq-v1.C \
                utility/src/sequence/dnaSeqFile-v1.C \
                utility/src/sequence/sequence-v1.C \
                \
                utility/src/system/logging-v1.C \
                utility/src/system/runtime-v1.C \
                utility/src/system/speedCounter-v1.C \
                utility/src/system/sweatShop-v1.C \
                utility/src/system/system-stackTrace-v1.C \
                utility/src/system/system-v1.C \
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
SOURCES      += utility/src/system/libbacktrace/atomic.c \
                utility/src/system/libbacktrace/backtrace.c \
                utility/src/system/libbacktrace/dwarf.c \
                utility/src/system/libbacktrace/elf.c \
                utility/src/system/libbacktrace/fileline.c \
                utility/src/system/libbacktrace/mmap.c \
                utility/src/system/libbacktrace/mmapio.c \
                utility/src/system/libbacktrace/posix.c \
                utility/src/system/libbacktrace/print.c \
                utility/src/system/libbacktrace/simple.c \
                utility/src/system/libbacktrace/sort.c \
                utility/src/system/libbacktrace/state.c \
                utility/src/system/libbacktrace/unknown.c
endif



SRC_INCDIRS  := . \
                utility/src \
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
                utgcns/utgcns.mk \
                utgcns/layoutToPackage.mk \
                \
                gfa/alignGFA.mk
