
.. _command-reference:

Canu Command Reference
======================

Every command, even the useless ones.

Commands marked as 'just usage' were automagically generated from the command line usage summary.  Yes, some of them even crashed.

:doc:`commands/bogart` (just usage)
  The unitig construction algorithm.  BOG stands for Best Overlap Graph; we haven't figured out what ART stands for.
:doc:`commands/bogus` (just usage)
  A unitig construction algorithm simulator.  Given reads mapped to a reference, returns the largest unitigs possible.
:doc:`commands/canu` (just usage)
  The executive in charge!  Coordinates all these commands to make an assembler.
:doc:`commands/correctOverlaps` (just usage)
  Part of Overlap Error Adjustment, recomputes overlaps given a set of read corrections.
:doc:`commands/estimate-mer-threshold` (just usage)
  Decides on a k-mer threshold for overlapInCore seeds.
:doc:`commands/fastqAnalyze` (just usage)
  Analyzes a FASTQ file and reports the best guess of the QV encoding.  Can also rewrite the FASTQ to be in Sanger QV format.
:doc:`commands/fastqSample` (just usage)
  Extracts random reads from a single or mated FASTQ file.  Extracts based on desired coverage, desired number of reads/pairs, desired fraction of the total, or desired total length.
:doc:`commands/fastqSimulate` (just usage)
  Creates reads with unbiased errors from a FASTA sequence.
:doc:`commands/fastqSimulate-sort` (just usage)
  Given input from fastqSimulate, sorts the reads by position in the reference.
:doc:`commands/filterCorrectionOverlaps` (just usage)
  Part of Read Correction, filters overlaps that shouldn't be used for correcting reads.
:doc:`commands/findErrors` (just usage)
  Part of Overlap Error Adjustment, generates a multialignment for each read, outputs a list of suspected errors in the read.
:doc:`commands/gatekeeperCreate` (just usage)
  Loads FASTA or FASTQ reads into the canu read database, gkpStore.
:doc:`commands/gatekeeperDumpFASTQ` (just usage)
  Outputs FASTQ reads fromt the canu read database, gkpStore.
:doc:`commands/gatekeeperDumpMetaData` (just usage)
  Outputs read and library metadata fromt the canu read database, gkpStore.
:doc:`commands/gatekeeperPartition` (just usage)
  Part of Consensus, rearranges the canu read database, gkpStore, to localize read to unitigs.
:doc:`commands/generateCorrectionLayouts` (just usage)
  Part of Read Correction, generates the multialignment layout used to correct reads.
:doc:`commands/leaff` (just usage)
  Not actually part of canu, but it came along with meryl.  Provides random access to FASTA, FASTQ and gkpStore.  Also does some analysis tasks.  Handy Swiss Army knife type of tool.
:doc:`commands/meryl` (just usage)
  Counts k-mer occurrences in FASTA, FASTQ and gkpStore.  Performs mathematical and logical operations on the resulting k-mer databases.
:doc:`commands/mhapConvert` (just usage)
  Convert mhap output to overlap output.
:doc:`commands/ovStoreBucketizer` (just usage)
  Part of the parallel overlap store building pipeline, loads raw overlaps from overlapper into the store.
:doc:`commands/ovStoreBuild` (just usage)
  Sequentially builds an overlap store from raw overlaps.  Simplest to run, but slow on large datasets.
:doc:`commands/ovStoreDump` (just usage)
  Dumps overlaps from the overlap store, ovlStore.
:doc:`commands/ovStoreIndexer` (just usage)
  Part of the parallel overlap store building pipeline, finalizes the store, after sorting with ovStoreSorter.
:doc:`commands/ovStoreSorter` (just usage)
  Part of the parallel overlap store building pipeline, sorts overlaps loaded into the store by ovStoreBucketizer.
:doc:`commands/overlapConvert` (just usage)
  Reads raw overlapper output, writes overlaps as ASCII.  The reverse of overlapImport.
:doc:`commands/overlapImport` (just usage)
  Reads ASCII overlaps in a few different formats, writes either 'raw overlapper output' or creates an ovlStore.
:doc:`commands/overlapInCore` (just usage)
  The classic overlapper algorithm.
:doc:`commands/overlapInCorePartition` (just usage)
  Generate partitioning to run overlapInCore jobs in parallel.
:doc:`commands/overlapPair` (just usage)
  An *experimental* algorithm to recompute overlaps and output the alignments.
:doc:`commands/prefixEditDistance-matchLimitGenerate` (just usage)
  Generate source code files with data representing the minimum length of a good overlap given some number of errors.
:doc:`commands/splitReads` (just usage)
  Part of Overlap Based Trimming, splits reads based on overlaps, specifically, looking for PacBio hairpin adapter signatures.
:doc:`commands/tgStoreCoverageStat` (just usage)
  Analyzes tigs in the tigStore, computes the classic `arrival rate statistic <http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Theory>`_.
:doc:`commands/tgStoreDump` (just usage)
  Analyzes and outputs tigs from the tigStore, in various formats (FASTQ, layouts, multialignments, etc).
:doc:`commands/tgStoreFilter` (just usage)
  Analyzes tigs in the tigStore, marks those that appear to be spurious 'degenerate' tigs.
:doc:`commands/tgStoreLoad` (just usage)
  Loads tigs into a tigStore.
:doc:`commands/tgTigDisplay` (just usage)
  Displays the tig contained in a binary multialignment file, as output by utgcns.
:doc:`commands/trimReads` (just usage)
  Part of Overlap Based Trimming, trims reads based on overlaps.
:doc:`commands/utgcns` (just usage)
  Generates a multialignment for a tig, based on the layout stored in tigStore.  Outputs FASTQ, layouts and binary mutlialignment files.
