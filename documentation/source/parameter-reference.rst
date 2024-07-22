
.. _parameter-reference:

Canu Parameter Reference
========================

To get the most up-to-date options, run

   canu -options

The default values below will vary based on the input data type and genome size.

Boolean options accept true/false or 1/0.

Memory sizes are assumed to be in gigabytes if no units are supplied.  Values may be non-integer
with or without a unit - 'k' for kilobytes, 'm' for megabytes, 'g' for gigabytes or 't' for
terabytes.  For example, "0.25t" is equivalent to "256g" (or simply "256").

Global Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The catch all category.

.. _rawErrorRate:

rawErrorRate <float=unset>
  The allowed difference in an overlap between two uncorrected reads, expressed as fraction error.
  Sets :ref:`corOvlErrorRate <corOvlErrorRate>` and :ref:`corErrorRate <corErrorRate>`.  The
  :ref:`rawErrorRate <rawErrorRate>` typically does not need to be modified.  It might need to be
  increased if very early reads are being assembled.  The default is 0.300 For PacBio reads, and
  0.500 for Nanopore reads.

.. _correctedErrorRate:

correctedErrorRate <float=unset>
  The allowed difference in an overlap between two corrected reads, expressed as fraction error.
  Sets :ref:`obtOvlErrorRate <obtOvlErrorRate>`, :ref:`utgOvlErrorRate <utgOvlErrorRate>`,
  :ref:`obtErrorRate <obtErrorRate>`, :ref:`utgErrorRate <utgErrorRate>`, and :ref:`cnsErrorRate
  <cnsErrorRate>`.
  The :ref:`correctedErrorRate <correctedErrorRate>` can be adjusted to account for the quality of
  read correction, for the amount of divergence in the sample being assembled, and for the amount of
  sequence being assembled.  The default is 0.045 for PacBio reads, and 0.144 for Nanopore reads.

  For low coverage datasets (less than 30X), we recommend increasing :ref:`correctedErrorRate
  <correctedErrorRate>` slightly, by 1% or so.

  For high-coverage datasets (more than 60X), we recommend decreasing :ref:`correctedErrorRate
  <correctedErrorRate>` slightly, by 1% or so.

  Raising the :ref:`correctedErrorRate <correctedErrorRate>` will increase run time.  Likewise,
  decreasing :ref:`correctedErrorRate <correctedErrorRate>` will decrease run time, at the risk of
  missing overlaps and fracturing the assembly.

.. _minReadLength:

minReadLength <integer=1000>
  Reads shorter than this are not loaded into the assembler.  Reads output by correction and
  trimming that are shorter than this are discarded.

  Must be no smaller than minOverlapLength.

  [In Canu v1.9 and earlier] If set high enough, the gatekeeper module will claim there are errors in the input reads,
  as too many of the input reads have been discarded.  As long as there is sufficient coverage,
  this is not a problem.  See :ref:`stopOnReadQuality <stopOnReadQuality>` and 
  :ref:`stopOnLowCoverage <stopOnLowCoverage>`

.. _minOverlapLength:

minOverlapLength <integer=500>
  Overlaps shorter than this will not be discovered.  Smaller values can be used to overcome lack of
  read coverage, but will also lead to false overlaps and potential misassemblies.  Larger values
  will result in more correct assemblies, but more fragmented, assemblies.

  Must be no bigger than minReadLength.

.. _readSamplingCoverage:

readSamplingCoverage <integer=unset>
  After loading all reads into the sequence store, discard some reads so that
  this amount of coverage remains.  Reads are discarded according to the score
  described in :ref:`readSamplingBias <readSamplingBias>`.

.. _readSamplingBias:

readSamplingBias <float=0.0>
  Adjust the sampling bias towards discarding longer (negative numbers) or
  shorter (positive numbers) reads.  Reads are assigned a score equal to
  `random_number * read_length ^ bias` and the lowest scoring reads are
  discarded, as described in :ref:`readSamplingCoverage
  <readSamplingCoverage>`.

  In the pictures below, green reads are kept, while purple reads are
  discarded.  The reads are along the X axis, sorted by decreasing score.
  The Y axis is the length of each read.

  A bias of 0.0 will retain random reads:

  .. image:: bias=+00.png

  A negative bias will retain shorter reads:

  .. image:: bias=-01.png
  .. image:: bias=-10.png

  A positive bias will retain longer reads:

  .. image:: bias=+01.png
  .. image:: bias=+10.png

.. _genomeSize:

genomeSize <float=unset> *required*
  An estimate of the size of the genome.  Common suffices are allowed, for example, 3.7m or 2.8g.

  The genome size estimate is used to decide how many reads to correct (via the :ref:`corOutCoverage <corOutCoverage>`
  parameter) and how sensitive the mhap overlapper should be (via the :ref:`mhapSensitivity <mhapSensitivity>`
  parameter). It also impacts some logging, in particular, reports of NG50 sizes.

.. _fast:

fast <toggle>
   This option uses MHAP overlapping for all steps, not just correction, making assembly significantly faster. It can be used on any genome size but may produce less continuous assemblies on genomes larger than 1 Gbp. It is recommended for nanopore genomes smaller than 1 Gbp or metagenomes.
   
   The fast option will also optionally use `wtdbg <https://github.com/ruanjue/wtdbg2>`_ for unitigging if wtdbg is manually copied to the Canu binary folder. However, this is only tested with very small genomes and is **NOT** recommended.

.. _canuIteration:

canuIteration <internal parameter, do not use>
  Which parallel iteration is being attempted.

canuIterationMax <integer=2>
  How many attempts, in total, Canu will make to execute a stage.  A stage
  is, e.g., computing overlaps for trimming, or computing consensus
  sequences, and is composed of several (typically hundreds) of compute jobs.

  The default setting of 2 allows a stage to fail once by *force majeure*
  (node failure, memory exaustion caused by other processes, etc) without
  stopping the entire workflow.  If a stage fails a second time the workflow
  will terminate.

  Any number of compute jobs can fail per attempt.  For example, if seven
  jobs on the first attempt failed, exactly those seven jobs will be executed
  again in the second attempt.

.. _onSuccess:

onSuccess <string=unset>
  Execute the command supplied when Canu successfully completes an assembly.  The command will
  execute in the <assembly-directory> (the -d option to canu) and will be supplied with the name of
  the assembly (the -p option to canu) as its first and only parameter.

.. _onFailure:

onFailure <string=unset>
  Execute the command supplied when Canu terminates abnormally.  The command will execute in the
  <assembly-directory> (the -d option to canu) and will be supplied with the name of the assembly
  (the -p option to canu) as its first and only parameter.

  There are two exceptions when the command is not executed: if a 'spec' file cannot be read, or if
  canu tries to access an invalid parameter.  The former will be reported as a command line error,
  and canu will never start.  The latter should never occur except when developers are developing
  the software.


Process Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _showNext:

showNext <boolean=false>
  Report the first major command that would be run, but don't run it.  Processing to get to that
  command, for example, checking the output of the previous command or preparing inputs for the next
  command, is still performed.

.. _stopOnReadQuality:

stopOnReadQuality <string=false>
  [In Canu v1.9 and earlier] If set, Canu will stop with the following error if there are significantly fewer reads or bases
  loaded into the read store than what is in the input data.

  ::

   Gatekeeper detected potential problems in your input reads.

   Please review the logging in files:
     /assembly/godzilla/asm.gkpStore.BUILDING.err
     /assembly/godzilla/asm.gkpStore.BUILDING/errorLog

   If you wish to proceed, rename the store with the following command and restart canu.

     mv /assembly/godzilla/asm.gkpStore.BUILDING \
        /assembly/godzilla/asm.gkpStore.ACCEPTED

   Option stopOnReadQuality=false skips these checks.

  The missing reads could be too short (decrease :ref:`minReadLength <minReadLength>` to include
  them), or have invalid bases or quality values.  A summary of the files loaded and errors detected
  is in the ``asm.gkpStore.BUILDING.err`` file, with full gory details in the
  ``asm.gkpStore.BUILDING/errorLog``.

  To proceed, set ``stopOnReadQuality=false`` or rename the directory as shown.

  Note that `U` bases are silently translated to `T` bases, to allow assembly of RNA sequences.

.. _stopOnLowCoverage:

stopOnLowCoverage <integer=10>
  Stop the assembly if read coverage is too low to be useful.  Coverage is
  checked whene when input sequences are
  initially loaded into the sequence store, when corrected reads are generated,
  and when read ends are trimmed off.


.. _stopAfter:

stopAfter <string=undefined>
  If set, Canu will stop processing after a specific stage in the pipeline finishes.  Valid values are:

  +-----------------------+-------------------------------------------------------------------+
  | **stopAfter=**        | **Canu will stop after ....**                                     |
  +-----------------------+-------------------------------------------------------------------+
  | sequenceStore         | reads are loaded into the assembler read database.                |
  +-----------------------+-------------------------------------------------------------------+
  | meryl-configure       | kmer counting jobs are configured.                                |
  +-----------------------+-------------------------------------------------------------------+
  | meryl-count           | kmers are counted, but not processed into one database.           |
  +-----------------------+-------------------------------------------------------------------+
  | meryl-merge           | kmers are merged into one database.                               |
  +-----------------------+-------------------------------------------------------------------+
  | meryl-process         | frequent kmers are generated.                                     |
  +-----------------------+-------------------------------------------------------------------+
  | meryl-subtract        | haplotype specific kmers are generated.                           |
  +-----------------------+-------------------------------------------------------------------+
  | meryl                 | all kmer work is complete.                                        |
  +-----------------------+-------------------------------------------------------------------+
  | haplotype-configure   | haplotype read separation jobs are configured.                    |
  +-----------------------+-------------------------------------------------------------------+
  | haplotype             | haplotype-specific reads are generated.                           |
  +-----------------------+-------------------------------------------------------------------+
  | overlapConfigure      | overlap jobs are configured.                                      |
  +-----------------------+-------------------------------------------------------------------+
  | overlap               | overlaps are generated, before they are loaded into the database. |
  +-----------------------+-------------------------------------------------------------------+
  | overlapStoreConfigure | the jobs for creating the overlap database are configured.        |
  +-----------------------+-------------------------------------------------------------------+
  | overlapStore          | overlaps are loaded into the overlap database.                    |
  +-----------------------+-------------------------------------------------------------------+
  | correction            | corrected reads are generated.                                    |
  +-----------------------+-------------------------------------------------------------------+
  | trimming              | trimmed reads are generated.                                      |
  +-----------------------+-------------------------------------------------------------------+
  | unitig                | unitigs and contigs are created.                                  |
  +-----------------------+-------------------------------------------------------------------+
  | consensusConfigure    | consensus jobs are configured.                                    |
  +-----------------------+-------------------------------------------------------------------+
  | consensus             | consensus sequences are loaded into the databases.                |
  +-----------------------+-------------------------------------------------------------------+

  *readCorrection* and *readTrimming* are deprecated synonyms for *correction* and *trimming*, respectively.

General Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

shell <string="/bin/sh">
  A path to a Bourne shell, to be used for executing scripts.  By default, '/bin/sh', which is typically
  the same as 'bash'.  C shells (csh, tcsh) are not supported.

java <string="java">
  A path to a Java application launcher of at least version 1.8.

minimap <string="minimap2">
  A path to the minimap2 versatile pairwise aligner.

gnuplot <string="gnuplot">
  A path to the gnuplot graphing utility.  Plotting is disabled if this is unset
  (`gnuplot=` or `gnuplot=undef`), or if gnuplot fails to execute, or if gnuplot
  cannot generate plots.  The latter two conditions generate warnings in the
  diagnostic output of Canu.

gnuplotImageFormat <string="png">
  The type of image to generate in gnuplot.  By default, canu will use png,
  svg or gif, in that order.

preExec <string=undef>
  A single command that will be run before Canu starts in a grid-enabled configuration.
  Can be used to set up the environment, e.g., with 'module'.

File Staging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The correction stage of Canu requires random access to all the reads.  Performance is greatly
improved if the gkpStore database of reads is copied locally to each node that computes corrected
read consensus sequences.  This 'staging' is enabled by supplying a path name to fast local storage
with the :ref:`stageDirectory <stageDirectory>` option, and, optionally, requesting access to that resource from the grid
with the :ref:`gridEngineStageOption <gridEngineStageOption>` option.

.. _stageDirectory:

stageDirectory <string=undefined>
  A path to a directory local to each compute node.  The directory should use an environment
  variable specific to the grid engine to ensure that it is unique to each task.

  For example, in Sun Grid Engine, `/scratch/$JOB_ID-$SGE_TASK_ID` will use both the numeric
  job ID and the numeric task ID.  In SLURM, `/scratch/$SLRUM_JOBID` accomplishes the same.

  If specified on the command line, be sure to escape the dollar sign, otherwise the shell will try
  to expand it before Canu sees the option: `stageDirectory=/scratch/\$JOB_ID-\$SGE_TASK_ID`.

  If specified in a specFile, do not escape the dollar signs.

.. _gridEngineStageOption:

gridEngineStageOption <string=undefined>
  This string is passed to the job submission command, and is expected to request
  local disk space on each node.  It is highly grid specific.  The string `DISK_SPACE`
  will be replaced with the amount of disk space needed, in gigabytes.

  On SLURM, an example is `--gres=lscratch:DISK_SPACE`


Cleanup Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

saveOverlaps <boolean=false>
  If 'true', retain all overlap stores.  If 'false', delete the correction
  and trimming overlap stores when they are no longer useful.  Overlaps used
  for contig construction are never deleted.

purgeOverlaps <string=normal>
  Controls when to remove intermediate overlap results.

  'never' removes no intermediate overlap results.  This is only useful if
  you have a desire to exhaust your disk space.

  'false' is the same as 'never'.

  'normal' removes intermediate overlap results after they are loaded into an
  overlap store.

  'true' is the same as 'normal'.

  'aggressive' removes intermediate overlap results as soon as possible.  In
  the event of a corrupt or lost file, this can result in a fair amount of
  suffering to recompute the data.  In particular, overlapper output is removed
  as soon as it is loaded into buckets, and buckets are removed once they are
  rewritten as sorted overlaps.

  'dangerous' removes intermediate results as soon as possible, in some
  cases, before they are even fully processed.  In addition to corrupt files,
  jobs killed by out of memory, power outages, stray cosmic rays, et cetera,
  will result in a fair amount of suffering to recompute the lost data.  This
  mode can help when creating ginormous overlap stores, by removing the
  bucketized data immediately after it is loaded into the sorting jobs, thus
  making space for the output of the sorting jobs.

  Use 'normal' for non-large assemblies, and when disk space is plentiful.
  Use 'aggressive' on large assemblies when disk space is tight.  Never use
  'dangerous', unless you know how to recover from an error and you fully
  trust your compute environment.

  For Mhap and Minimap2, the raw ovelraps (in Mhap and PAF format) are
  deleted immediately after being converted to Canu ovb format, except when
  purgeOverlaps=never.

saveReadCorrections <boolean=false>.
  If set, do not remove raw corrected read output from correction/2-correction. Normally, this
  output is removed once the corrected reads are generated.

saveIntermediates <boolean=false>
  If set, do not remove intermediate outputs.  Normally, intermediate files are removed
  once they are no longer needed.

  NOT IMPLEMENTED.

saveMerCounts <boolean=false>
  If set, do not remove meryl binary databases.

saveReads <boolean=false>
  If set, save the corrected reads (in asm.correctedReads.fasta.gz) and trimmed reads (in asm.trimmedReads.fasta.gz).
  Both read sets are saved in the asm.gkpStore, and can be retrieved later.

Executive Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Canu 'executive' is responsible for controlling what tasks run and when they run.  It doesn't
directly do any significant computations, rather it just examines the files that exist and decides
which component to run next.  For example, if overlaps exist but contigs do not, it would create
contigs next.

When under grid control, some tasks can be run in the same job as the executive, if there is enough
memory and threads reserved for the executive.  The benefit of this is slight; on a heavily loaded
grid, it would reduce the number of job scheduling iterations Canu needs to run.

.. _executiveMemory:

executiveMemory <integer=4>

  The amount of memory, in gigabytes, to reserve when running the Canu exectuve (and any jobs it
  runs directly).  Increasing this past 4 GB can allow some tasks (such as creating an overlap store
  or creating contigs) to run directly, without needing a separate grid job.

.. _executiveThreads:

executiveThreads <integer=1>

  The number of threads to reserve for the Canu executive.


Overlapper Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Overlaps are generated for three purposes: read correction, read trimming and unitig construction.
The algorithm and parameters used can be set independently for each set of overlaps.

Two overlap algorithms are in use.  One, mhap, is typically applied to raw uncorrected reads and
returns alignment-free overlaps with imprecise extents.  The other, the original overlapper
algorithm 'ovl', returns alignments but is much more expensive.

There are three sets of parameters, one for the 'mhap' algorithm, one for the 'ovl' algorithm, and
one for the 'minimap' algorithm.  Parameters used for a specific type of overlap are set by a prefix
on the option: 'cor' for read correction, 'obt' for read trimming ('overlap based trimming') or
'utg' for unitig construction.  For example, 'corOverlapper=ovl' would set the overlapper used for
read correction to the 'ovl' algorithm.

{prefix}Overlapper <string=see-below>
  Specify which overlap algorithm, 'mhap' or 'ovl' or 'minimap'.  The default is to use 'mhap' for
  'cor' and 'ovl' for both 'obt' and 'utg'.

Overlapper Configuration, ovl Algorithm
---------------------------------------

.. _corOvlErrorRate:
.. _obtOvlErrorRate:
.. _utgOvlErrorRate:
.. _ovlErrorRate:

{prefix}OvlErrorRate <float=unset>
  Overlaps above this error rate are not computed.
  * :ref:`corOvlErrorRate <corOvlErrorRate>` applies to overlaps generated for correcting reads;
  * :ref:`obtOvlErrorRate <obtOvlErrorRate>` applied to overlaps generated for trimming reads;
  * :ref:`utgOvlErrorRate <utgOvlErrorRate>` applies to overlaps generated for assembling reads.
  These limits apply to the 'ovl' overlap algorithm and when alignments are computed for mhap
  overlaps with :ref:`reAlign <reAlign>`.

{prefix}OvlFrequentMers <string=undefined>
  Do not seed overlaps with these kmers, or, for mhap, do not seed with these kmers unless necessary (down-weight them).

  For corFrequentMers (mhap), the file must contain a single line header followed by number-of-kmers data lines::

    0 number-of-kmers
    forward-kmer word-frequency kmer-count total-number-of-kmers
    reverse-kmer word-frequency kmer-count total-number-of-kmers

  Where `kmer-count` is the number of times this kmer sequence occurs in the reads, 'total-number-of-kmers'
  is the number of kmers in the reads (including duplicates; rougly the number of bases in the reads),
  and 'word-frequency' is 'kmer-count' / 'total-number-of-kmers'.

  For example::

    0 4
    AAAATAATAGACTTATCGAGTC  0.0000382200    52      1360545
    GACTCGATAAGTCTATTATTTT  0.0000382200    52      1360545
    AAATAATAGACTTATCGAGTCA  0.0000382200    52      1360545
    TGACTCGATAAGTCTATTATTT  0.0000382200    52      1360545

  This file must be gzip compressed.

  For obtFrequentMers and ovlFrequentMers, the file must contain a list of the canonical kmers and
  their count on a single line.  The count value is ignored, but needs to be present.  This file
  should not be compressed.

  For example::

    AAAATAATAGACTTATCGAGTC  52
    AAATAATAGACTTATCGAGTCA  52

{prefix}OvlHashBits <integer=unset>
  Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per ovlHashBlockLength.

{prefix}OvlHashBlockLength <integer=unset>
  Amount of sequence (bp to load into the overlap hash table.

{prefix}OvlHashLoad <integer=unset>
  Maximum hash table load.  If set too high, table lookups are inefficient; if too low, search
  overhead dominates run time.

{prefix}OvlMerDistinct <integer=unset>
  K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps.

{prefix}OvlMerSize <integer=unset>
  K-mer size for seeds in overlaps.

{prefix}OvlMerThreshold <integer=unset>
  K-mer frequency threshold; mers more frequent than this count are not used to seed overlaps.

{prefix}OvlMerTotal <integer=unset>
  K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps.

{prefix}OvlRefBlockLength <integer=unset>
  Amount of sequence (bp to search against the hash table per batch.

{prefix}OvlRefBlockSize <integer=unset>
  Number of reads to search against the hash table per batch.

Overlapper Configuration, mhap Algorithm
----------------------------------------

{prefix}MhapBlockSize <integer=unset>
  For the MHAP overlapper, the number of reads to load per GB of memory (mhapMemory).
  When mhapSensitivity=high, this value is automatically divided by two.

{prefix}MhapMerSize <integer=unset>
  K-mer size for seeds in mhap.

.. _mhapReAlign:

{prefix}ReAlign <boolean=false>
  Compute actual alignments from mhap overlaps.
  uses either obtErrorRate or ovlErrorRate, depending on which overlaps are computed)

.. _mhapSensitivity:

{prefix}MhapSensitivity <string="normal">
  Coarse sensitivity level: 'low', 'normal' or 'high'.  Based on read coverage (which is impacted by
  genomeSize), 'low' sensitivity is used if coverage is more than 60; 'normal' is used if coverage
  is between 60 and 30, and 'high' is used for coverages less than 30.

Overlapper Configuration, mmap Algorithm
----------------------------------------

{prefix}MMapBlockSize <integer=unset>
  Number of reads per 1GB block.  Memory * size is loaded into memory per job.

{prefix}MMapMerSize <integer=unset>
  K-mer size for seeds in minimap.

Overlap Store
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The overlap algorithms return overlaps in an arbitrary order, however, all other algorithms (or
nearly all) require all overlaps for a single read to be readily available.  Thus, the overlap store
collects and sorts the overlapper outputs into a store of overlaps, sorted by the first read in the
overlap.  Each overlap is listed twice in the store, once in an "A vs B" format, and once in a "B vs
A" format (that is, swapping which read is 'first' in the overlap description).

Two construction algorithms are supported.  A 'sequential' method uses a single data stream, and is
faster for small and moderate size assemblies.  A 'parallel' method uses multiple compute nodes and
can be faster (depending on your network disk bandwidth) for moderate and large assemblies.  Be
advised that the parallel method is less efficient than the sequential method, and can easily thrash
consumer-level NAS devices resulting in exceptionally poor performance.

The sequential method load all overlapper outputs (.ovb files in 1-overlapper) into memory,
duplicating each overlap.  It then sortes overlaps, and creates the final overlap store.

The parallel method uses two parallel tasks: bucketizing ('ovb' tasks) and sorting ('ovs' tasks).
Bucketizing reads the outputs of the overlap tasks (ovb files in 1-overlapper), duplicates each
overlap, and writes these to intermediate files.  Sorting tasks load these intermediate file into
memory, sorts the overlaps, then writes the sorted overlaps back to disk.  There will be one
'bucketizer' ('ovb' tasks) task per overlap task, and tens to hundreds of 'sorter' ('ovs' tasks).  A
final 'indexing' step is done in the Canu executive, which ties all the various files togather into
the final overlap store.

Increasing ovsMemory will allow more overlaps to fit into memory at once.  This will allow larger
assemblies to use the sequential method, or reduce the number of 'ovs' tasks for the parallel
method.

Increasing the allowed memory for the Canu executive can allow the overlap store to be constructed as
part of the executive job -- a separate grid job for constructing the store is not needed.

ovsMemory <float>
  How much memory, in gigabytes, to use for constructing overlap stores.  Must be at least 256m or 0.25g.

Meryl
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 'meryl' algorithm counts the occurrences of kmers in the input reads.  It outputs a FASTA format
list of frequent kmers, and (optionally) a binary database of the counts for each kmer in the input.

Meryl can run in (almost) any memory size, by splitting the computation into smaller (or larger) chunks.

merylMemory <integer=unset>
  Amount of memory, in gigabytes, to use for counting kmers.

merylThreads <integer=unset>
  Number of compute threads to use for kmer counting.


Overlap Based Trimming
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _obtErrorRate:

obtErrorRate <float=unset>
  Stringency of overlaps to use for trimming reads.

trimReadsOverlap <integer=1>
  Minimum overlap between evidence to make contiguous trim.

trimReadsCoverage <integer=1>
  Minimum depth of evidence to retain bases.


Trio binning Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _hapUnknownFraction:

hapUnknownFraction <float=0.05>
  Fraction of unclassified bases to ignore for haplotype assemblies. If there are more than this fraction of unclassified bases, they are included in both haplotype assemblies.

.. _grid-engine:

Grid Engine Support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Canu directly supports most common grid scheduling systems.  Under normal use, Canu will query the
system for grid support, configure itself for the machines available in the grid, then submit itself
to the grid for execution.  The Canu pipeline is a series of about a dozen steps that alternate
between embarrassingly parallel computations (e.g., overlap computation) and sequential bookkeeping
steps (e.g., checking if all overlap jobs finished).  This is entirely managed by Canu.

Canu has first class support for the various schedulers derived from Sun Grid Engine (Univa, Son of
Grid Engine) and the Simple Linux Utility for Resource Management (SLURM), meaning that the
developers have direct access to these systems.  Platform Computing's Load Sharing Facility (LSF) and
the various schedulers derived from the Portable Batch System (PBS, Torque and PBSPro) are supported
as well, but without developer access bugs do creep in.  As of Canu v1.5, support seems stable and
working.

useGrid <boolean=true>
  Master control.  If 'false', no algorithms will run under grid control.  Does not change the value
  of the other useGrid options.

  If 'remote', jobs are configured for grid execution, but not submitted.  A message, with commands
  to launch the job, is reported and canu halts execution.

  Note that the host used to run canu for 'remote' execution must know about the grid, that is, it
  must be able to submit jobs to the grid.

It is also possible to enable/disable grid support for individual algorithms with options such as
`useGridBAT`, `useGridCNS`, et cetera.  This has been useful in the (far) past to prevent certain
algorithms, notably overlap error adjustment, from running too many jobs concurrently and thrashing
disk.  Recent storage systems seem to be able to handle the load better -- computers have gotten
faster quicker than genomes have gotten larger.

There are many options for configuring a new grid ('gridEngine*') and for configuring how canu
configures its computes to run under grid control ('gridOptions*').  The grid engine to use is
specified with the 'gridEngine' option.

gridEngine <string>
  Which grid engine to use.  Auto-detected.  Possible choices are 'sge', 'pbs', 'pbspro', 'lsf' or
  'slurm'.

.. _grid-engine-config:

Grid Engine Configuration
-------------------------

There are many options to configure support for a new grid engine, and we don't describe them fully.
If you feel the need to add support for a new engine, please contact us.  That said, file
``src/pipeline/canu/Defaults.pm`` lists a whole slew of parameters that are used to build up grid
commands, they all start with ``gridEngine``.  For each grid, these parameters are defined in the
various ``src/pipeline/Grid_*.pm`` modules.  The parameters are used in
``src/pipeline/canu/Execution.pm``.

In Canu 1.8 and earlier, ``gridEngineMemoryOption`` and ``gridEngineThreadsOption`` are used to tell
Canu how to request resources from the grid.  Starting with ``snapshot v1.8 +90 changes`` (roughly
January 11th), those options were merged into ``gridEngineResourceOption``.  These options specify
the grid options needed to request memory and threads for each job.  For example, the default
``gridEngineResourceOption`` for PBS/Torque is "-l nodes=1:ppn=THREADS:mem=MEMORY", and for Slurm it
is "--cpus-per-task=THREADS --mem-per-cpu=MEMORY".  Canu will replace "THREADS" and "MEMORY" with
the specific values needed for each job.

.. _grid-options:

Grid Options
------------

To run on the grid, each stage needs to be configured - to tell the grid how many cores it will use and how much memory it needs.
Some support for this is automagic (for example, overlapInCore and mhap know how to do this), others need to be manually configured.
Yes, it's a problem, and yes, we want to fix it.

The gridOptions* parameters supply grid-specific options to the grid submission command.

gridOptions <string=unset>
  Grid submission command options applied to all grid jobs
gridOptionsJobName <string=unset>
  Grid submission command jobs name suffix
gridOptionsBAT <string=unset>
  Grid submission command options applied to unitig construction with the bogart algorithm
gridOptionsGFA <string=unset>
  Grid submission command options applied to gfa alignment and processing
gridOptionsCNS <string=unset>
  Grid submission command options applied to unitig consensus jobs
gridOptionsCOR <string=unset>
  Grid submission command options applied to read correction jobs
gridOptionsExecutive <string=unset>
  Grid submission command options applied to master script jobs
gridOptionsOEA <string=unset>
  Grid submission command options applied to overlap error adjustment jobs
gridOptionsRED <string=unset>
  Grid submission command options applied to read error detection jobs
gridOptionsOVB <string=unset>
  Grid submission command options applied to overlap store bucketizing jobs
gridOptionsOVS <string=unset>
  Grid submission command options applied to overlap store sorting jobs
gridOptionsCORMHAP <string=unset>
  Grid submission command options applied to mhap overlaps for correction jobs
gridOptionsCOROVL <string=unset>
  Grid submission command options applied to overlaps for correction jobs
gridOptionsOBTMHAP <string=unset>
  Grid submission command options applied to mhap overlaps for trimming jobs
gridOptionsOBTOVL <string=unset>
  Grid submission command options applied to overlaps for trimming jobs
gridOptionsUTGMHAP <string=unset>
  Grid submission command options applied to mhap overlaps for unitig construction jobs
gridOptionsUTGOVL <string=unset>
  Grid submission command options applied to overlaps for unitig construction jobs



Algorithm Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several algorithmic components of canu can be disabled, based on the type of the reads being
assembled, the type of processing desired, or the amount of compute resources available.  Overlap

enableOEA <boolean=true>
  Do overlap error adjustment - comprises two steps: read error detection (RED and overlap error adjustment (OEA

Algorithm Execution Method
--------------------------

Canu has a fairly sophisticated (or complicated, depending on if it is working or not) method for
dividing large computes, such as read overlapping and consensus, into many smaller pieces and then
running those pieces on a grid or in parallel on the local machine.  The size of each piece is
generally determined by the amount of memory the task is allowed to use, and this memory size --
actually a range of memory sizes -- is set based on the genomeSize parameter, but can be set
explicitly by the user.  The same holds for the number of processors each task can use.
For example, a genomeSize=5m would result in overlaps using between 4gb and
8gb of memory, and between 1 and 8 processors.

Given these requirements, Canu will pick a specific memory size and number of processors
so that the maximum number of jobs will run at the same time.  In the overlapper example,
if we are running on a machine with 32gb memory and 8 processors, it is not possible to run
8 concurrent jobs that each require 8gb memory, but it is possible to run 4 concurrent jobs
each using 6gb memory and 2 processors.

To completely specify how Canu runs algorithms, one needs to specify a maximum memory size, a
maximum number of processors, and how many pieces to run at one time.  Users can set these manually
through the {prefix}Memory, {prefix}Threads and {prefix}Concurrency options.  If they are not
set, defaults are chosen based on genomeSize.

{prefix}Concurrency <integer=unset>
  Set the number of tasks that can run at the same time, when running without grid support.

{prefix}Threads <integer=unset>
  Set the number of compute threads used per task.

{prefix}Memory <integer=unset>
  Set the amount of memory, in gigabytes, to use for each job in a task.

Available prefixes are:

+-------+-----------+---------------------------------------------------------------+
|    Prefix         | Algorithm                                                     |
+=======+===========+===============================================================+
| cor   |           | | Overlap generation using the 'mhap' algorithm for           |
+-------+           | | 'cor'=correction                                            |
| obt   | mhap      | | 'obt'=trimming                                              |
+-------+           | | 'utg'=assembly                                              |
| utg   |           |                                                               |
+-------+-----------+---------------------------------------------------------------+
| cor   |           | | Overlap generation using the 'minimap2' algorithm for       |
+-------+           | | 'cor' = correction                                          |
| obt   | mmap      | | 'obt' = trimming                                            |
+-------+           | | 'utg' = assembly                                            |
| utg   |           |                                                               |
+-------+-----------+---------------------------------------------------------------+
| cor   |           | | Overlap generation using the 'overlapInCore' algorithm for  |
+-------+           | | 'cor' = correction                                          |
| obt   | ovl       | | 'obt' = trimming                                            |
+-------+           | | 'utg' = assembly                                            |
| utg   |           |                                                               |
+-------+-----------+---------------------------------------------------------------+
|       | ovb       | Parallel overlap store bucketizing                            |
+-------+-----------+---------------------------------------------------------------+
|       | ovs       | Parallel overlap store bucket sorting                         |
+-------+-----------+---------------------------------------------------------------+
|       | cor       | Read correction                                               |
+-------+-----------+---------------------------------------------------------------+
|       | red       | Error detection in reads                                      |
+-------+-----------+---------------------------------------------------------------+
|       | oea       | Error adjustment in overlaps                                  |
+-------+-----------+---------------------------------------------------------------+
|       | bat       | Unitig/contig construction                                    |
+-------+-----------+---------------------------------------------------------------+
|       | cns       | Unitig/contig consensus                                       |
+-------+-----------+---------------------------------------------------------------+

For example, 'mhapMemory` would set the memory limit for computing overlaps with the mhap algorithm;
'cormhapMemory' would set the memory limit only when mhap is used for generating overlaps used for
correction.

.. _maxMemory:
.. _minMemory:
.. _minThreads:
.. _maxThreads:

The 'minMemory', 'maxMemory', 'minThreads' and 'maxThreads' options will apply to all jobs, and
can be used to artificially limit canu to a portion of the current machine.  In the overlapper
example above, setting maxThreads=4 would result in two concurrent jobs instead of four.


Overlap Error Adjustment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Overlap Error Adjustment module adjusts the error rate claimed by each
overlap to account for sequencing error and true polymorphism.  A pair-wise
multialignment is generated for all overlaps to a given read.  Each
multialignment column is examined to determine if the base in the given read
is correct, is part of a true difference, or is a likely sequencing error.
For the latter case, a base change is noted.  Once all base changes in all
columns of all reads have been determined, all overlaps are recomputed with
said changes applied, and the new fraction error stored for each overlap.

Three parameters exist to change the behavior:

oeaErrorRate <real=unset>
  Do not use overlaps that are above this fraction error for finding true
  differences and errors.  Mostly useful only when oeaMaskTrivial is enabled.
  For PacBio and Nanopore data , the default is utgOvlErrorRate.  For HiFi
  data, the default is 0.003.

oeaHaploConfirm <integer=5>
  A set of this many reads with the same base in a column is declared to be a
  true haplotype difference and will increase the fraction error for the
  overlap.  The default for all read types is 5.
  
oeaMaskTrivial <boolean=unset>
  If true, errors that occur in low-complexity sequence are ignored when
  finding base changes and when recomputing overlap error rates.  Unless set
  explicitly, it is true (enabled) for PacBio HiFi data and false otherwise.
  
This module consists of two steps: RED (read error detection), and OEA
(overlap error adjustment).  They have slightly different run-time
requirements.  RED can use multiple threads and is slightly more
computationally expensive; OEA can not use multiple threads and is slightly
more I/O intensive.  The batch length and batch size parameters can tune the
size of each job, however, the default values have worked well (so well that
we don't actually document what these values should be set to).

oeaBatchLength <unset>
  Number of bases per overlap error correction batch
oeaBatchSize <unset>
  Number of reads per overlap error correction batch

redBatchLength <unset>
  Number of bases per fragment error detection batch
redBatchSize <unset>
  Number of reads per fragment error detection batch


Unitigger
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unitigger <string="bogart">
  Which unitig construction algorithm to use.  Only "bogart" is supported.

.. _utgErrorRate:

utgErrorRate <float=unset>
  Stringency of overlaps used for constructing contigs.  The `bogart` algorithm uses the distribution of overlap error rates to filter high error overlaps; `bogart` will never see overlaps with error higher than this parameter.

batOptions <unset>
  Advanced options to bogart

Consensus Partitioning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

STILL DONE BY UNITIGGER, NEED TO MOVE OUTSIDE

cnsConsensus
  Which algorithm to use for computing consensus sequences.  Only 'utgcns' is supported.

cnsPartitions
  Compute conseus by splitting the tigs into N partitions.

cnsPartitionMin
  Don't make a partition with fewer than N reads

cnsMaxCoverage
  Limit unitig consensus to at most this coverage.

.. _cnsErrorRate:

cnsErrorRate
  Inform the consensus generation algorithm of the amount of difference it should expect in a
  read-to-read alignment.  Typically set to :ref:`utgOvlErrorRate <utgOvlErrorRate>`.  If set too
  high, reads could be placed in an incorrect location, leading to errors in the consensus sequence.
  If set too low, reads could be omitted from the consensus graph (or multialignment, depending on
  algorithm), resulting in truncated consensus sequences.

.. _correction:

Read Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in Canu is to find high-error overlaps and generate corrected sequences for
subsequent assembly. This is currently the fastest step in Canu. By default, only the longest 40X of
data (based on the specified genome size) is used for correction. Typically, some reads are trimmed
during correction due to being chimeric or having erroneous sequence, resulting in a loss of 20-25%
(30X output). You can force correction to be non-lossy by setting `corMinCoverage=0`, in which case
the corrected reads output will be the same length as the input data, keeping any high-error
unsupported bases. Canu will trim these in downstream steps before assembly.

If you have a dataset with uneven coverage or small plasmids, correcting the longest 40X may not
give you sufficient coverage of your genome/plasmid. In these cases, you can set
`corOutCoverage=999`, or any value greater than your total input coverage which will correct and
assemble all input data, at the expense of runtime.

.. _corErrorRate:

corErrorRate <integer=unset>
  Do not use overlaps with error rate higher than this when computing corrected reads.

  In Canu v2.2, this parameter was changed from 0.50 (for -nanopore) and 0.30
  (for -pacbio) to 0.30 and 0.25, respectively.

  The tables below show a significant speedup for Nanopore reads without much
  loss in output quantity.  There is indication of a slight improvement in
  corrected read quality at lower corErrorRate, however, read quality was not
  directly evaluated.

  For PacBio reads, with a smaller change in corErrorRate, the speedup is
  about 10%.

  +-----------------+--------------------+------------------+-------------------+------------+
  | CHM13 Chromosome X, nanopore, 105x input coverage                                        |
  +-----------------+--------------------+------------------+-------------------+------------+
  | |               | | Corrected        | | Trimmed        | | Bogart          | | CPU Time |
  | | corErrorRate  | | Coverage         | | Coverage       | | Error Rate      | | (hours)  |
  +=================+====================+==================+===================+============+
  | 5               | 22.0x              | 21.8x            | 0.3958%           |            |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 10              | 35.7x              | 35.2x            |                   |            |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 15              | 38.1x              | 37.5x            |                   |            |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 20              | 38.6x              | 38.0x            |                   | 1160       |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 25              | 38.7x              | 38.1x            |                   | 1290       |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 30              | 38.8x              | 38.1x            |                   | 1449       |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 40              | 38.8x              | 38.1x            |                   | 1625       |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 50              | 38.8x              | 38.2x            |                   | 3683       |
  +-----------------+--------------------+------------------+-------------------+------------+

  +-----------------+--------------------+------------------+-------------------+------------+
  | HG002 Chromosome X, nanopore, 20x input coverage                                         |
  +-----------------+--------------------+------------------+-------------------+------------+
  | |               | | Corrected        | | Trimmed        | | Bogart          | | CPU Time |
  | | corErrorRate  | | Coverage         | | Coverage       | | Error Rate      | | (hours)  |
  +=================+====================+==================+===================+============+
  | 5               | --.-x              | --.-x            | -%                | 31         |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 10              |  3.9x              | --.-x            | -%                | 66         |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 15              |  9.6x              | --.-x            | -%                | 105        |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 20              | 11.4x              | 11.2x            | 1.71%             | 134        |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 25              | 11.9x              | 11.6x            | 1.79%             | 154        |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 30              | 12.0x              | 11.8x            | 1.83%             | 169        |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 40              | 12.2x              | 12.0x            | 1.94%             | 221        |
  +-----------------+--------------------+------------------+-------------------+------------+
  | 50              | 12.6x              | 12.3x            | 2.29%             | 709        |
  +-----------------+--------------------+------------------+-------------------+------------+

corConsensus <string="falconpipe">
  Which algorithm to use for computing read consensus sequences.  Only 'falcon' and 'falconpipe' are supported.

corPartitions <integer=128>
  Partition read correction into N jobs

corPartitionMin <integer=25000>
  Don't make a read correction partition with fewer than N reads

corMinEvidenceLength <integer=unset>
  Limit read correction to only overlaps longer than this; default: unlimited
corMinCoverage <integer=4>
  Limit read correction to regions with at least this minimum coverage. Split reads when coverage drops below threshold.
corMaxEvidenceErate <integer=unset>
  Limit read correction to only overlaps at or below this fraction error; default: unlimited
corMaxEvidenceCoverageGlobal <string="1.0x">
  Limit reads used for correction to supporting at most this coverage; default: 1.0 * estimated coverage
corMaxEvidenceCoverageLocal <string="2.0x">
  Limit reads being corrected to at most this much evidence coverage; default: 10 * estimated coverage

.. _corOutCoverage:

corOutCoverage <integer=40>
  Only correct the longest reads up to this coverage; default 40

corFilter <string="expensive">
  Method to filter short reads from correction; 'quick' or 'expensive' or 'none'

Output Filtering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _contigFilter:

contigFilter <minReads, integer=2> <minLength, integer=0> <singleReadSpan, float=1.0> <lowCovSpan, float=0.5> <lowCovDepth, integer=5>
  A contig that needs any of the following conditions is flagged as 'unassembled' and removed from
  further consideration:

    - fewer than minReads reads (default 2)
    - shorter than minLength bases (default 0)
    - a single read covers more than singleReadSpan fraction of the contig (default 1.0)
    - more than lowCovSpan fraction of the contig is at coverage below lowCovDepth (defaults 0.5, 5)

  This filtering is done immediately after initial contigs are formed, before potentially
  incorrectly spanned repeats are detected.  Initial contigs that incorrectly span a repeat can be
  split into multiple contigs; none of these new contigs will be flagged as 'unassembled', even if
  they are a single read.
