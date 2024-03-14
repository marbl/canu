MODULE       :=    canu
VERSION      := snapshot v2.2
VERSION_H    := utility/src/version.H

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
                utility/src/parasail/cigar.c \
                utility/src/parasail/memory.c \
                utility/src/parasail/sg.c \
                utility/src/parasail/sg_trace.c \
                utility/src/parasail/sg_qb_de_dispatch.c \
                utility/src/parasail/sg_qe_db_dispatch.c \
                utility/src/parasail/sg_qx_dispatch.c \
                utility/src/parasail/sg_trace.c \
                \
                utility/src/sequence/dnaSeq-v1.C \
                utility/src/sequence/bufSeqFile-v1.C \
                utility/src/sequence/htsSeqFile-v1.C \
                utility/src/sequence/sequence-v1.C \
                \
                utility/src/system/logging-v1.C \
                utility/src/system/runtime-v1.C \
                utility/src/system/speedCounter-v1.C \
                utility/src/system/sweatShop-v1.C \
                utility/src/system/system-stackTrace-v1.C \
                utility/src/system/system-v1.C \
                utility/src/system/time-v1.C \
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
                gfa/gfa.C \
                gfa/bed.C

SOURCES      += utility/src/htslib/hts/bcf_sr_sort.c \
                utility/src/htslib/hts/bgzf.c \
                utility/src/htslib/hts/errmod.c \
                utility/src/htslib/hts/faidx.c \
                utility/src/htslib/hts/header.c \
                utility/src/htslib/hts/hfile.c \
                utility/src/htslib/hts/hfile_libcurl.c \
                utility/src/htslib/hts/hfile_s3.c \
                utility/src/htslib/hts/hts.c \
                utility/src/htslib/hts/hts_expr.c \
                utility/src/htslib/hts/hts_os.c \
                utility/src/htslib/hts/kfunc.c \
                utility/src/htslib/hts/kstring.c \
                utility/src/htslib/hts/md5.c \
                utility/src/htslib/hts/multipart.c \
                utility/src/htslib/hts/probaln.c \
                utility/src/htslib/hts/realn.c \
                utility/src/htslib/hts/regidx.c \
                utility/src/htslib/hts/region.c \
                utility/src/htslib/hts/sam.c \
                utility/src/htslib/hts/synced_bcf_reader.c \
                utility/src/htslib/hts/tbx.c \
                utility/src/htslib/hts/textutils.c \
                utility/src/htslib/hts/thread_pool.c \
                utility/src/htslib/hts/vcf.c \
                utility/src/htslib/hts/vcf_sweep.c \
                utility/src/htslib/hts/vcfutils.c \
                utility/src/htslib/cram/cram_codecs.c \
                utility/src/htslib/cram/cram_decode.c \
                utility/src/htslib/cram/cram_encode.c \
                utility/src/htslib/cram/cram_external.c \
                utility/src/htslib/cram/cram_index.c \
                utility/src/htslib/cram/cram_io.c \
                utility/src/htslib/cram/cram_stats.c \
                utility/src/htslib/cram/mFILE.c \
                utility/src/htslib/cram/open_trace_file.c \
                utility/src/htslib/cram/pooled_alloc.c \
                utility/src/htslib/cram/string_alloc.c \
                utility/src/htslib/htscodecs/arith_dynamic.c \
                utility/src/htslib/htscodecs/fqzcomp_qual.c \
                utility/src/htslib/htscodecs/htscodecs.c \
                utility/src/htslib/htscodecs/pack.c \
                utility/src/htslib/htscodecs/rANS_static.c \
                utility/src/htslib/htscodecs/rANS_static32x16pr.c \
                utility/src/htslib/htscodecs/rANS_static32x16pr_neon.c \
                utility/src/htslib/htscodecs/rANS_static4x16pr.c \
                utility/src/htslib/htscodecs/rle.c \
                utility/src/htslib/htscodecs/tokenise_name3.c \
                utility/src/htslib/htscodecs/utils.c


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

SRC_INCDIRS  := utility/src \
                stores \
                stores/libsnappy \
                overlapInCore \
                overlapInCore/liboverlap

SYS_INCDIRS  += $(shell pkg-config --cflags-only-I openssl libcurl liblzma | sed s:-I/:/:g)
LDFLAGS      += $(shell pkg-config --libs-only-L   openssl libcurl liblzma)
LDLIBS       += $(shell pkg-config --libs-only-l   openssl libcurl liblzma) -lz -lbz2

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
                seqrequester/src/seqrequester/seqrequester.mk \
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

EXECUTABLES := pipelines/canu.pl                        -> bin/canu \
               pipelines/canu-time.pl                   -> bin/canu-time \
               pipelines/draw-tig.pl                    -> bin/draw-tig

FILES       += pipelines/canu.defaults                  -> bin/canu.defaults \
               mhap/mhap-2.1.3.jar                      -> share/java/classes/mhap-2.1.3.jar \
               pipelines/canu/Consensus.pm              -> lib/perl5/site_perl/canu/Consensus.pm \
               pipelines/canu/CorrectReads.pm           -> lib/perl5/site_perl/canu/CorrectReads.pm \
               pipelines/canu/HaplotypeReads.pm         -> lib/perl5/site_perl/canu/HaplotypeReads.pm \
               pipelines/canu/Configure.pm              -> lib/perl5/site_perl/canu/Configure.pm \
               pipelines/canu/Defaults.pm               -> lib/perl5/site_perl/canu/Defaults.pm \
               pipelines/canu/Execution.pm              -> lib/perl5/site_perl/canu/Execution.pm \
               pipelines/canu/SequenceStore.pm          -> lib/perl5/site_perl/canu/SequenceStore.pm \
               pipelines/canu/Grid.pm                   -> lib/perl5/site_perl/canu/Grid.pm \
               pipelines/canu/Grid_Cloud.pm             -> lib/perl5/site_perl/canu/Grid_Cloud.pm \
               pipelines/canu/Grid_DNANexus.pm          -> lib/perl5/site_perl/canu/Grid_DNANexus.pm \
               pipelines/canu/Grid_LSF.pm               -> lib/perl5/site_perl/canu/Grid_LSF.pm \
               pipelines/canu/Grid_Local.pm             -> lib/perl5/site_perl/canu/Grid_Local.pm \
               pipelines/canu/Grid_PBSTorque.pm         -> lib/perl5/site_perl/canu/Grid_PBSTorque.pm \
               pipelines/canu/Grid_SGE.pm               -> lib/perl5/site_perl/canu/Grid_SGE.pm \
               pipelines/canu/Grid_Slurm.pm             -> lib/perl5/site_perl/canu/Grid_Slurm.pm \
               pipelines/canu/Meryl.pm                  -> lib/perl5/site_perl/canu/Meryl.pm \
               pipelines/canu/Output.pm                 -> lib/perl5/site_perl/canu/Output.pm \
               pipelines/canu/OverlapBasedTrimming.pm   -> lib/perl5/site_perl/canu/OverlapBasedTrimming.pm \
               pipelines/canu/OverlapErrorAdjustment.pm -> lib/perl5/site_perl/canu/OverlapErrorAdjustment.pm \
               pipelines/canu/OverlapInCore.pm          -> lib/perl5/site_perl/canu/OverlapInCore.pm \
               pipelines/canu/OverlapMhap.pm            -> lib/perl5/site_perl/canu/OverlapMhap.pm \
               pipelines/canu/OverlapMMap.pm            -> lib/perl5/site_perl/canu/OverlapMMap.pm \
               pipelines/canu/OverlapStore.pm           -> lib/perl5/site_perl/canu/OverlapStore.pm \
               pipelines/canu/Report.pm                 -> lib/perl5/site_perl/canu/Report.pm \
               pipelines/canu/Unitig.pm                 -> lib/perl5/site_perl/canu/Unitig.pm
