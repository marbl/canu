

.. _celera-assembler: `Celera Assembler <http://wgs-assembler.sourceforge.net>`

.. _tutorial:

Canu Tutorial
=============

Canu assembles reads from PacBio RS II or Oxford Nanopore MinION instruments into
uniquely-assemblable contigs, unitigs.  Canu owes lots of it design and code to
`celera-assembler <Celera Assembler>`_.

Canu can be run using hardware of nearly any shape or size, anywhere from laptops to computational
grids with thousands of nodes.  Obviously, larger assemblies will take a long time to compute on
laptops, and smaller assemblies can't take advantage of hundreds of nodes, so what is being
assembled plays some part in determining what hardware can be effectively used.

Most algorithms in canu have been multi-threaded (to use all the cores on a single node),
parallelized (to use all the nodes in a grid), or both (all the cores on all the nodes).

.. _canu-command:

Canu, the command
~~~~~~~~~~~~~~~~~~~~~~

The **canu** command is the 'executive' program that runs all modules of the assembler.  It oversees
each of the three top-level tasks (correction, trimming, unitig construction), each of which
consists of many steps.  Canu ensures that input files for each step exist, that each step
successfully finished, and that the output for each step exists.  It does minor bits of processing,
such as reformatting files, but generally just executes other programs.

::

 canu [-trimmed|-untrimmed|-raw|-corrected] \
   [-s <assembly-specifications-file>] \
    -p <assembly-prefix> \
    -d <assembly-directory> \
    genomeSize=<number>[g|m|k] \
    [other-options] \
    [-pacbio|-nanopore|-pacbio-hifi] *fastq

The -p option, to set the file name prefix of intermediate and output files, is mandatory.  If -d is
not supplied, canu will run in the current directory, otherwise, Canu will create the
`assembly-directory` and run in that directory.  It is _not_ possible to run two different
assemblies in the same directory.

The -s option will import a list of parameters from the supplied specification ('spec') file.  These
parameters will be applied before any from the command line are used, providing a method for
setting commonly used parameters, but overriding them for specific assemblies.

By default, all needed top-level tasks are performed (-pacbio and -nanopore are assumed to be raw and untrimmed while -pacbio-hifi are assumed to be corrected and trimmed).  It is possible to run exactly one task by
specifying your read characteristics.  These options can be useful if you want to correct
reads once and try many different assemblies.  We do exactly that in the :ref:`quickstart`.
Additionally, suppling pre-corrected reads with -pacbio -corrected or -nanopore -corrected
will run only the trimming and assembling stages. Specifying reads as -corrected -untrimmed will run only the assembly step.

Parameters are key=value pairs that configure the assembler.  They set run time parameters (e.g.,
memory, threads, grid), algorithmic parameters (e.g., error rates, trimming aggressiveness), and
enable or disable entire processing steps (e.g., don't correct errors, don't search for subreads).
They are described later.  One parameter is required: the genomeSize (in bases, with common
SI prefixes allowed, for example, 4.7m or 2.8g; see :ref:`genomeSize`).  Parameters are listed in
the :ref:`parameter-reference`, but the common ones are described in this document.

Reads are supplied to canu by options that options that describe how the reads were generated, and
what level of quality they are, for example, -pacbio-raw indicates the reads were generated on a
PacBio RS II instrument, and have had no processing done to them.  Each file of reads supplied this
way becomes a 'library' of reads.  The reads should have been (physically) generated all at the same
time using the same steps, but perhaps sequenced in multiple batches.  In canu, each library has a
set of options setting various algorithmic parameters, for example, how aggressively to trim.  To
explicitly set library parameters, a text 'gkp' file describing the library and the input files must
be created.  Don't worry too much about this yet, it's an advanced feature, fully described in
Section :ref:`gkp-files`.

The read-files contain sequence data in either FASTA or FASTQ format (or both!  A quirk of the
implementation allows files that contain both FASTA and FASTQ format reads).  The files can be
uncompressed, gzip, bzip2 or xz compressed.  We've found that "gzip -1" provides good compression
that is fast to both compress and decompress.  For 'archival' purposes, we use "xz -9".

.. _canu-pipeline:

Canu, the pipeline
~~~~~~~~~~~~~~~~~~~~~~

The canu pipeline, that is, what it actually computes, comprises of computing overlaps and
processing the overlaps to some result.  Each of the three tasks (read correction, read trimming and
unitig construction) follow the same pattern:

* Load reads into the read database, gkpStore.
* Compute k-mer counts in preparation for the overlap computation.
* Compute overlaps.
* Load overlaps into the overlap database, ovlStore.
* Do something interesting with the reads and overlaps.

  * The read correction task will replace the original noisy read sequences with consensus sequences
    computed from overlapping reads.
  * The read trimming task will use overlapping reads to decide what regions of each read are
    high-quality sequence, and what regions should be trimmed.  After trimming, the single largest
    high-quality chunk of sequence is retained.
  * The unitig construction task finds sets of overlaps that are consistent, and uses those to place
    reads into a multialignment layout.  The layout is then used to generate a consensus sequence
    for the unitig.



.. _module-tags:

Module Tags
~~~~~~~~~~~~~~~~~~~~~~

Because each of the three tasks share common algorithms (all compute overlaps, two compute
consensus sequences, etc), parameters are differentiated by a short prefix 'tag' string.  This lets
canu have one generic parameter that can be set to different values for each stage in each task.
For example, "corOvlMemory" will set memory usage for overlaps being generated for read correction;
"obtOvlMemory" for overlaps generated for Overlap Based Trimming; "utgOvlMemory" for overlaps
generated for unitig construction.

The tags are:

+--------+-------------------------------------------------------------------+
|Tag     | Usage                                                             |
+========+===================================================================+
|master  | the canu script itself, and any components that it runs directly  |
+--------+-------------------------------------------------------------------+
+--------+-------------------------------------------------------------------+
|cns     | unitig consensus generation                                       |
+--------+-------------------------------------------------------------------+
|cor     | read correction generation                                        |
+--------+-------------------------------------------------------------------+
+--------+-------------------------------------------------------------------+
|red     | read error detection                                              |
+--------+-------------------------------------------------------------------+
|oea     | overlap error adjustment                                          |
+--------+-------------------------------------------------------------------+
+--------+-------------------------------------------------------------------+
|ovl     | the standard overlapper                                           |
+--------+-------------------------------------------------------------------+
|corovl  | the standard overlapper, as used in the correction phase          |
+--------+-------------------------------------------------------------------+
|obtovl  | the standard overlapper, as used in the trimming phase            |
+--------+-------------------------------------------------------------------+
|utgovl  | the standard overlapper, as used in the assembly phase            |
+--------+-------------------------------------------------------------------+
+--------+-------------------------------------------------------------------+
|mhap    | the mhap overlapper                                               |
+--------+-------------------------------------------------------------------+
|cormhap | the mhap overlapper, as used in the correction phase              |
+--------+-------------------------------------------------------------------+
|obtmhap | the mhap overlapper, as used in the trimming phase                |
+--------+-------------------------------------------------------------------+
|utgmhap | the mhap overlapper, as used in the assembly phase                |
+--------+-------------------------------------------------------------------+
+--------+-------------------------------------------------------------------+
|mmap    | the `minimap <https://github.com/lh3/minimap>`_ overlapper        |
+--------+-------------------------------------------------------------------+
|cormmap | the minimap overlapper, as used in the correction phase           |
+--------+-------------------------------------------------------------------+
|obtmmap | the minimap overlapper, as used in the trimming phase             |
+--------+-------------------------------------------------------------------+
|utgmmap | the minimap overlapper, as used in the assembly phase             |
+--------+-------------------------------------------------------------------+
+--------+-------------------------------------------------------------------+
|ovb     | the bucketizing phase of overlap store building                   |
+--------+-------------------------------------------------------------------+
|ovs     | the sort phase of overlap store building                          |
+--------+-------------------------------------------------------------------+

We'll get to the details eventually.

.. _execution:

Execution Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

There are two modes that canu runs in: locally, using just one machine, or grid-enabled, using
multiple hosts managed by a grid engine.  LSF, PBS/Torque, PBSPro, Sun Grid Engine (and
derivations), and Slurm are supported, though LSF has had limited testing. Section
:ref:`grid-engine-config` has a few hints on how to set up a new grid engine.

By default, if a grid is detected the canu pipeline will immediately submit itself to the grid and
run entirely under grid control.  If no grid is detected, or if option ``useGrid=false`` is set,
canu will run on the local machine.

In both cases, Canu will auto-detect available resources and configure job sizes based on the
resources and genome size you're assembling. Thus, most users should be able to run the command
without modifying the defaults. Some advanced options are outlined below. Each stage has the same
five configuration options, and tags are used to specialize the option to a specific stage.  The
options are:

useGrid<tag>=boolean
  Run this stage on the grid, usually in parallel.
gridOptions<tag>=string
  Supply this string to the grid submit command.
<tag>Memory=integer
  Use this many gigabytes of memory, per process.
<tag>Threads
  Use this many compute threads per process.
<tag>Concurrency
  If not on the grid, run this many jobs at the same time.

Global grid options, applied to every job submitted to the grid, can be set with 'gridOptions'.
This can be used to add accounting information or access credentials.

A name can be associated with this compute using 'gridOptionsJobName'.  Canu will work just fine
with no name set, but if multiple canu assemblies are running at the same time, they will tend to
wait for each others jobs to finish.  For example, if two assemblies are running, at some point both
will have overlap jobs running.  Each assembly will be waiting for all jobs named 'ovl_asm' to
finish.  Had the assemblies specified job names, gridOptionsJobName=apple and
gridOptionsJobName=orange, then one would be waiting for jobs named 'ovl_asm_apple', and the other
would be waiting for jobs named 'ovl_asm_orange'.

.. _error-rates:

Error Rates
~~~~~~~~~~~~~~~~~~~~~~

Canu expects all error rates to be reported as fraction error, not as percent error.  We're not sure
exactly why this is so.  Previously, it used a mix of fraction error and percent error (or both!),
and was a little confusing.  Here's a handy table you can print out that converts between fraction
error and percent error.  Not all values are shown (it'd be quite a large table) but we have every
confidence you can figure out the missing values:

==============  =============
Fraction Error  Percent Error
==============  =============
0.01            1%
0.02            2%
0.03            3%
.               .
.               .
0.12            12%
.               .
.               .
==============  =============

Canu error rates always refer to the percent difference in an alignment of two reads, not the
percent error in a single read, and not the amount of variation in your reads.  These error rates
are used in two different ways: they are used to limit what overlaps are generated, e.g., don't
compute overlaps that have more than 5% difference; and they are used to tell algorithms what
overlaps to use, e.g., even though overlaps were computed to 5% difference, don't trust any above 3%
difference.

There are seven error rates.  Three error rates control overlap creation (:ref:`corOvlErrorRate
<ovlErrorRate>`, :ref:`obtOvlErrorRate <ovlErrorRate>` and :ref:`utgOvlErrorRate <ovlErrorRate>`),
and four error rates control algorithms (:ref:`corErrorRate <corErrorRate>`, :ref:`obtErrorRate
<obtErrorRate>`, :ref:`utgErrorRate <utgErrorRate>`, :ref:`cnsErrorRate <cnsErrorRate>`).

The three error rates for overlap creation apply to the `ovl` overlap algorithm and the
:ref:`mhapReAlign <mhapReAlign>` option used to generate alignments from `mhap` or `minimap`
overlaps.  Since `mhap` is used for generating correction overlaps, the :ref:`corOvlErrorRate
<ovlErrorRate>` parameter is not used by default.  Overlaps for trimming and assembling use the
`ovl` algorithm, therefore, :ref:`obtOvlErrorRate <ovlErrorRate>` and :ref:`utgOvlErrorRate
<ovlErrorRate>` are used.

The four algoriothm error rates are used to select which overlaps can be used for correcting reads
(:ref:`corErrorRate <corErrorRate>`); which overlaps can be used for trimming reads
(:ref:`obtErrorRate <obtErrorRate>`); which overlaps can be used for assembling reads
(:ref:`utgErrorRate <utgErrorRate>`).  The last error rate, :ref:`cnsErrorRate <cnsErrorRate>`,
tells the consensus algorithm to not trust read alignments above that value.

For convenience, two meta options set the error rates used with uncorrected reads
(:ref:`rawErrorRate <rawErrorRate>`) or used with corrected reads.  (:ref:`correctedErrorRate
<correctedErrorRate>`).  The default depends on the type of read being assembled.

================== ======  ========
Parameter          PacBio  Nanopore
================== ======  ========
rawErrorRate       0.300   0.500
correctedErrorRate 0.045   0.144
================== ======  ========

In practice, only :ref:`correctedErrorRate <correctedErrorRate>` is usually changed.  The :ref:`faq`
has :ref:`specific suggestions <tweak>` on when to change this.

Canu v1.4 and earlier used the :ref:`errorRate <errorRate>` parameter, which set the expected
rate of error in a single corrected read.

.. _minimum-lengths:

Minimum Lengths
~~~~~~~~~~~~~~~~~~~~~~

Two minimum sizes are known:

minReadLength
  Discard reads shorter than this when loading into the assembler, and when trimming reads.

minOverlapLength
  Do not save overlaps shorter than this.

Overlap configuration
~~~~~~~~~~~~~~~~~~~~~~

The largest compute of the assembler is also the most complicated to configure.  As shown in the
'module tags' section, there are up to eight (!) different overlapper configurations.  For
each overlapper ('ovl' or 'mhap') there is a global configuration, and three specializations
that apply to each stage in the pipeline (correction, trimming or assembly).

Like with 'grid configuration', overlap configuration uses a 'tag' prefix applied to each option.  The
tags in this instance are 'cor', 'obt' and 'utg'.

For example:

- To change the k-mer size for all instances of the ovl overlapper, 'merSize=23' would be used.
- To change the k-mer size for just the ovl overlapper used during correction, 'corMerSize=16' would be used.
- To change the mhap k-mer size for all instances, 'mhapMerSize=18' would be used.
- To change the mhap k-mer size just during correction, 'corMhapMerSize=15' would be used.
- To use minimap for overlap computation just during correction, 'corOverlapper=minimap' would be used. The minimap2 executable must be symlinked from the Canu binary folder ('Linux-amd64/bin' or 'Darwin-amd64/bin' depending on your system).

Ovl Overlapper Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

<tag>Overlapper
  select the overlap algorithm to use, 'ovl' or 'mhap'.

Ovl Overlapper Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

<tag>ovlHashBlockLength
  how many bases to reads to include in the hash table; directly controls process size
<tag>ovlRefBlockSize
  how many reads to compute overlaps for in one process; directly controls process time
<tag>ovlRefBlockLength
 same, but use 'bases in reads' instead of 'number of reads'
<tag>ovlHashBits
  size of the hash table (SHOULD BE REMOVED AND COMPUTED, MAYBE TWO PASS)
<tag>ovlHashLoad
  how much to fill the hash table before computing overlaps (SHOULD BE REMOVED)
<tag>ovlMerSize
  size of kmer seed; smaller - more sensitive, but slower

The overlapper will not use frequent kmers to seed overlaps.  These are computed by the 'meryl' program,
and can be selected in one of three ways.

Terminology.  A k-mer is a contiguous sequence of k bases.  The read 'ACTTA' has two 4-mers: ACTT
and CTTA.  To account for reverse-complement sequence, a 'canonical kmer' is the lexicographically
smaller of the forward and reverse-complemented kmer sequence.  Kmer ACTT, with reverse complement
AAGT, has a canonical kmer AAGT.  Kmer CTTA, reverse-complement TAAG, has canonical kmer CTTA.

A 'distinct' kmer is the kmer sequence with no count associated with it.  A 'total' kmer (for lack
of a better term) is the kmer with its count.  The sequence TCGTTTTTTTCGTCG has 12 'total' 4-mers
and 8 'distinct' kmers.

::

 TCGTTTTTTTCGTCG  count
 TCGT             2 distinct-1
  CGTT            1 distinct-2
   GTTT           1 distinct-3
    TTTT          4 distinct-4
     TTTT         4 copy of distinct-4
      TTTT        4 copy of distinct-4
       TTTT       4 copy of distinct-4
        TTTC      1 distinct-5
         TTCG     1 distinct-6
          TCGT    2 copy of distinct-1
           CGTC   1 distinct-7
            GTCG  1 distinct-8


<tag>MerThreshold
  any kmer with count higher than N is not used
<tag>MerDistinct
  pick a threshold so as to seed overlaps using this fraction of all distinct kmers in the input.  In the example above,
  fraction 0.875 of the k-mers (7/8) will be at or below threshold 2.
<tag>MerTotal
  pick a threshold so as to seed overlaps using this fraction of all kmers in the input.  In the example above,
  fraction 0.667 of the k-mers (8/12) will be at or below threshold 2.
<tag>FrequentMers
  don't compute frequent kmers, use those listed in this file

Mhap Overlapper Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

<tag>MhapBlockSize
  Chunk of reads that can fit into 1GB of memory. Combined with memory to compute the size of chunk the reads are split into.
<tag>MhapMerSize
  Use k-mers of this size for detecting overlaps.
<tag>ReAlign
  After computing overlaps with mhap, compute a sequence alignment for each overlap.
<tag>MhapSensitivity
  Either 'normal', 'high', or 'fast'.

Mhap also will down-weight frequent kmers (using tf-idf), but it's selection of frequent is not exposed.

Minimap Overlapper Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~
<tag>MMapBlockSize
  Chunk of reads that can fit into 1GB of memory. Combined with memory to compute the size of chunk the reads are split into.
<tag>MMapMerSize
  Use k-mers of this size for detecting overlaps

Minimap also will ignore high-frequency minimizers, but it's selection of frequent is not exposed.

.. _outputs:

Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As Canu runs, it outputs status messages, execution logs, and some analysis to the console.  Most of
the analysis is captured in ``<prefix>.report`` as well.

LOGGING

<prefix>.report
  Most of the analysis reported during assembly. This will report the histogram of read lengths, the histogram or k-mers in the raw and corrected reads, the summary of corrected data, summary of overlaps, and the summary of contig lengths. 
  
  You can use the k-mer corrected read histograms with tools like `GenomeScope <http://qb.cshl.edu/genomescope/>`_ to estimate heterozygosity and genome size. In particular, histograms with more than 1 peak likely indicate a heterozygous genome. See the :ref:`FAQ` for some suggested parameters.
  
  The corrected read report gives a summary of the fate of all input reads. The first part:::

    --                             original      original
    --                            raw reads     raw reads
    --   category                w/overlaps  w/o/overlaps
    --   -------------------- ------------- -------------
    --   Number of Reads             250609           477
    --   Number of Bases         2238902045       1896925
    --   Coverage                    97.344         0.082
    --   Median                        6534          2360
    --   Mean                          8933          3976
    --   N50                          11291          5756
    --   Minimum                       1012             0
    --   Maximum                      60664         41278

  reports the fraction of reads which had an overlap. In this case, the majority had at least one overlap, which is good. Next::
  
    --                                        --------corrected---------  
    --                             evidence                     expected      
    --   category                     reads            raw     corrected
    --   -------------------- -------------  ------------- -------------
    --   Number of Reads             229397          48006         48006      
    --   Number of Bases         2134291652      993586222     920001699     
    --   Coverage                    92.795         43.199        40.000          
    --   Median                        6842          15330         14106           
    --   Mean                          9303          20697         19164           
    --   N50                          11512          28066         26840           
    --   Minimum                       1045          10184         10183         
    --   Maximum                      60664          60664         59063     
    --   

  reports that a total of 92.8x of raw bases are candidates for correction. By default, Canu only selects the longest 40x for correction. In this case, it selects 43.2x of raw read data which it estimates will result in 40x correction. Not all raw reads survive full-length through correction::

    --                            ----------rescued----------
    --                                        expected
    --   category                     raw     corrected
    --   --------------------   ------------- -------------
    --   Number of Reads               20030         20030
    --   Number of Bases            90137165      61903752
    --   Coverage                      3.919         2.691
    --   Median                         3324          2682
    --   Mean                           4500          3090
    --   N50                            5529          3659
    --   Minimum                        1012           501
    --   Maximum                       41475         10179

  The rescued reads are those which would not have contributed to the correction of the selected longest 40x subset. These could be short plasmids, mitochondria, etc. Canu includes them even though they're too short by the 40x cutoff to avoid losing sequence during assembly. Lastly::

    --                        --------uncorrected--------
    --                                           expected
    --   category                       raw     corrected
    --   -------------------- ------------- -------------
    --   Number of Reads             183050        183050
    --   Number of Bases         1157075583     951438105
    --   Coverage                    50.308        41.367
    --   Median                        5729          5086
    --   Mean                          6321          5197
    --   N50                           7467          6490
    --   Minimum                          0             0
    --   Maximum                      50522         10183

  are the reads which were deemed too short to correct. If you increase ``corOutCoverage``, you could get up to 41x more corrected sequence. However, unless the genome is very heterozygous, this does not typically improve the assembly and increases the running time.

  The assembly read error report summarizes how unitigging was run:::
  
    --  ERROR RATES
    --  -----------
    --                                                   --------threshold------
    --  21017                        fraction error      fraction        percent
    --  samples                              (1e-5)         error          error
    --                   --------------------------      --------       --------
    --  command line (-eg)                           ->     30.00        0.0300%
    --  command line (-eM)                           ->   1000.00        1.0000%
    --  mean + std.dev       0.08 +-   4 *     1.05  ->      4.30        0.0043%
    --  median + mad         0.00 +-   4 *     0.00  ->      0.00        0.0000%
    --  90th percentile                              ->      1.00        0.0010%  (enabled)
    --  

  Canu selects multiple error rate thresholds and selects the most appropriate one based on how many reads end up without overlaps at each threshold. In this case, it used 0.001% or 1 error in 10 kbp after considering 1.0% and 0.0043%.
  
  The assembly statistics (NG50, etc) are reported before and after consensus calling. Note that for HiFi data, the pre-consensus statistics are in homopolymer-compressed space.

READS

<prefix>.correctedReads.fasta.gz
   The reads after correction.

<prefix>.trimmedReads.fasta.gz
   The corrected reads after overlap based trimming.

SEQUENCE

<prefix>.contigs.fasta
   Everything which could be assembled and is the full assembly, including both unique, repetitive, and bubble elements.

<prefix>.unassembled.fasta
   Reads and low-coverage contigs which could not be incorporated into the primary assembly.

The header line for each sequence provides some metadata on the sequence.::

   >tig######## len=<integer> reads=<integer> covStat=<float> class=<contig|bubble|unassm> suggestRepeat=<yes|no> suggestCircular=<yes|no>

   len
      Length of the sequence, in bp.

   reads
      Number of reads used to form the contig.

   covStat
      The log of the ratio of the contig being unique versus being two-copy, based on the read arrival rate.  Positive values indicate more likely to be unique, while negative values indicate more likely to be repetitive.  See `Footnote 24 <http://science.sciencemag.org/content/287/5461/2196.full#ref-24>`_ in `Myers et al., A Whole-Genome Assembly of Drosophila <http://science.sciencemag.org/content/287/5461/2196.full>`_.

   class
      Type of sequence.  Unassembled sequences are primarily low-coverage sequences spanned by a single read.

   suggestRepeat
      If yes, sequence was detected as a repeat based on graph topology or read overlaps to other sequences.

   suggestCircular
      If yes, sequence is likely circular.  The fasta def line includes the non-redundant coordinates

GRAPHS

Canu versions prior to v1.9 created a GFA of the contig graph.  However, as noted at the time, the
GFA format cannot represent partial overlaps between contigs (for more details see the discussion of
general edges on the `GFA2 <https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md>`_ page).
Because Canu contigs are not compatible with the GFA format, <prefix>.contigs.gfa has been removed. Prior to Canu v2.1, 
contigs split at overlap junctions were output as unitigs. However, these graphs often would be missing edges and
be over-fragmented (split where there is no ambiguity). Thus <prefix>.unitigs.fasta and <prefix.unitigs.gfa have been removed.

METADATA

The layout provides information on where each read ended up in the final assembly, including
contig and positions. It also includes the consensus sequence for each contig.

<prefix>.contigs.layout
  (undocumented)

<prefix>.contigs.layout.readToTig
  The position of each read in a contig.
  
  The file looks like::
  
    #readID	tigID	coordType	bgn	end
    677083	4	ungapped	0	23436
    2343812	4	ungapped	12469	1223

  In this case read ids 677083 and 2343812 ended up in tig00000004 and the coordinates are listed at the end (read 2343812 is reverse-complemented). 
  
  You need to do a bit of work to get the original id of 2343812, look in the gkpStore/readNames.txt file, there you should find::

   2343812	m54033_180126_223601/39780749/39781_51526 id=4778961_0 id=2354708 clr=181,11399

  which gives you the original read (PacBio in this case) id.

<prefix>.contigs.layout.tigInfo
  A list of the contigs, lengths, coverage, number of reads and other metadata.
  Essentially the same information provided in the FASTA header line.






