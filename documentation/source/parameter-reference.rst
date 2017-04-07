
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

.. _errorRate:

errorRate <float=unset> (OBSOLETE)
  This parameter was removed on January 27th, 2016, and is valid only in Canu 1.4 or earlier.
  Canu currently still accepts the :ref:`errorRate <errorRate>` parameter, but its use is strongly discouraged.

  The expected error in a single corrected read.  The seven error rates were then set to three times
  this value (except for :ref:`corErrorRate <corErrorRate>`).

.. _rawErrorRate:

rawErrorRate <float=unset>
  The allowed difference in an overlap between two uncorrected reads, expressed as fraction error.
  Sets :ref:`corOvlErrorRate` and :ref:`corErrorRate`.  The `rawErrorRate` typically does not need
  to be modified.  It might need to be increased if very early reads are being assembled.  The
  default is 0.300 For PacBio reads, and 0.500 for Nanopore reads.

.. _correctedErrorRate:

correctedErrorRate <float=unset>
  The allowed difference in an overlap between two corrected reads, expressed as fraction error.  Sets :ref:`obtOvlErrorRate`, :ref:`utgOvlErrorRate`, :ref:`obtErrorRate`, :ref:`utgErrorRate`, and :ref:`cnsErrorRate`.
  The `correctedErrorRate` can be adjusted to account for the quality of read correction, for the amount of divergence in the sample being
  assembled, and for the amount of sequence being assembled.  The default is 0.045 for PacBio reads, and 0.144 for Nanopore reads.

  For low coverage datasets (less than 30X), we recommend increasing `correctedErrorRate` slightly, by 1% or so.

  For high-coverage datasets (more than 60X), we recommend decreasing `correctedErrorRate` slighly, by 1% or so.

  Raising the `correctedErrorRate` will increase run time.  Likewise, decreasing `correctedErrorRate` will decrease run time, at the risk of missing overlaps and fracturing the assembly.

.. _minReadLength:

minReadLength <integer=1000>
  Reads shorter than this are not loaded into the assembler.  Reads output by correction and
  trimming that are shorter than this are discarded.

  Must be no smaller than minOverlapLength.

  If set high enough, the gatekeeper module will halt as too many of the input reads have been
  discarded.  Set `stopOnReadQuality` to false to avoid this.

.. _minOverlapLength:

minOverlapLength <integer=500>
  Overlaps shorter than this will not be discovered.  Smaller values can be used to overcome lack of
  read coverage, but will also lead to false overlaps and potential misassemblies.  Larger values
  will result in more correct assemblies, but more fragmented, assemblies.

  Must be no bigger than minReadLength.

.. _genomeSize:

genomeSize <float=unset> *required*
  An estimate of the size of the genome.  Common suffices are allowed, for example, 3.7m or 2.8g.

  The genome size estimate is used to decide how many reads to correct (via the corOutCoverage_
  parameter) and how sensitive the mhap overlapper should be (via the mhapSensitivity_
  parameter). It also impacts some logging, in particular, reports of NG50 sizes.

.. _canuIteration:

canuIteration <internal parameter, do not use>
  Which parallel iteration is being attempted.
canuIterationMax <integer=2>
  How many parallel iterations to try.  Ideally, the parallel jobs, run under grid control, would all finish successfully on the first try.
  Sometimes, jobs fail due to other jobs exhausting resources (memory), or by the node itself failing.  In this case, canu will launch the jobs
  again.  This parameter controls how many times it tries.

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

showNext <boolean=false>
  Report the first major command that would be run, but don't run it.  Processing to get to that
  command, for example, checking the output of the previous command or preparing inputs for the
  next command, is still performed.

stopAfter <string=undefined>
  If set, Canu will stop processing after a specific stage in the pipeline finishes.

  Valid values for ``stopAfter`` are:

   - ``gatekeeper`` - stops after the reads are loaded into the assembler read database.
   - ``meryl`` - stops after frequent kmers are tabulated.
   - ``overlapConfigure`` - stops after overlap jobs are configured.
   - ``overlap`` - stops after overlaps are generated, before they are loaded into the overlap database.
   - ``overlapStoreConfigure`` - stops after the ``ovsMethod=parallel`` jobs are configured; has no impact for ``ovsMethod=sequential``.
   - ``overlapStore`` - stops after overlaps are loaded into the overlap database.
   - ``readCorrection`` - stops after corrected reads are generated.
   - ``readTrimming`` - stops after trimmed reads are generated.
   - ``unitig`` - stops after unitigs and contigs are created.
   - ``consensusConfigure`` - stops after consensus jobs are configured.
   - ``consensus`` - stops after consensus sequences are loaded into the databases.


General Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pathMap <string=undefined>
  A text file containing lines that map a hostname to a path to the assembler binaries.
  This can be used to provide fine-grained binary directories, for example, when incompatible versions
  of the operating system are used, or when canu is installed in a non-standard way.

  The hostname must be the same as returned by 'uname -n'.  For example::

   grid01   /grid/software/canu/Linux-amd64/bin/
   devel01  /devel/canu/Linux-amd64/bin/

shell <string="/bin/sh">
  A path to a Bourne shell, to be used for executing scripts.  By default, '/bin/sh', which is typically
  the same as 'bash'.  C shells (csh, tcsh) are not supported.

java <string="java">
  A path to a Java application launcher of at least version 1.8.

gnuplot <string="gnuplot">
  A path to the gnuplot graphing utility.

gnuplotImageFormat <string="png">
  The type of image to generate in gnuplot.  By default, canu will use png, svg or gif, in that order.

gnuplotTested <boolean=false>
  If set, skip the tests to determine if gnuplot will run, and to decide the image type to generate.  This is used when gnuplot fails to run, or isn't even installed, and allows canu to continue execution without generating graphs.


File Staging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The correction stage of Canu requires random access to all the reads.  Performance is greatly
improved if the gkpStore database of reads is copied locally to each node that computes corrected
read consensus sequences.  This 'staging' is enabled by supplying a path name to fast local storage
with the `stageDirectory` option, and, optionally, requesting access to that resource from the grid
with the `gridEngineStageOption` option.

stageDirectory <string=undefined>
  A path to a directory local to each compute node.  The directory should use an environment
  variable specific to the grid engine to ensure that it is unique to each task.

  For example, in Sun Grid Engine, `/scratch/$JOB_ID-$SGE_TASK_ID` will use both the numeric
  job ID and the numeric task ID.  In SLURM, `/scratch/$SLRUM_JOBID` accomplishes the same.

  If specified on the command line, be sure to escape the dollar sign, otherwise the shell will try
  to expand it before Canu sees the option: `stageDirectory=/scratch/\$JOB_ID-\$SGE_TASK_ID`.

  If specified in a specFile, do not escape the dollar signs.

gridEngineStageOption <string=undefined>
  This string is passed to the job submission command, and is expected to request
  local disk space on each node.  It is highly grid specific.  The string `DISK_SPACE`
  will be replaced with the amount of disk space needed, in gigabytes.

  On SLURM, an example is `--gres=lscratch:DISK_SPACE`


Cleanup Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

saveOverlaps <boolean=false>
  If set, do not remove raw overlap output from either mhap or overlapInCore.  Normally, this output is removed once
  the overlaps are loaded into an overlap store.

saveReadCorrections <boolean=false.
  If set, do not remove raw corrected read output from correction/2-correction. Normally, this output is removed once the corrected reads are generated.
  
saveIntermediates <boolean=false>
  If set, do not remove intermediate outputs.  Normally, intermediate files are removed
  once they are no longer needed.

  NOT IMPLEMENTED.

saveMerCounts <boolean=false>
  If set, do not remove meryl binary databases.

Overlapper Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Overlaps are generated for three purposes: read correction, read trimming and unitig construction.
The algorithm and parameters used can be set independently for each set of overlaps.

Two overlap algorithms are in use.  One, mhap, is typically applied to raw uncorrected reads and
returns alignment-free overlaps with imprecise extents.  The other, the original overlapper
algorithm 'ovl', returns alignments but is much more expensive.

There are three sets of parameters, one for the 'mhap' algorithm, one for the 'ovl' algorithm, and one for the 'minimap' algorithm.
Parameters used for a specific type of overlap are set by a prefix on the option: 'cor' for read
correction, 'obt' for read trimming ('overlap based trimming') or 'utg' for unitig construction.
For example, 'corOverlapper=ovl' would set the overlapper used for read correction to the 'ovl'
algorithm.

{prefix}Overlapper <string=see-below>
  Specify which overlap algorith, 'mhap' or 'ovl' or 'minimap'.  The default is to use 'mhap' for 'cor' and 'ovl' for both 'obt' and 'utg'.

Overlapper Configuration, ovl Algorithm
---------------------------------------

.. _corOvlErrorRate:
.. _obtOvlErrorRate:
.. _utgOvlErrorRate:
.. _ovlErrorRate:

{prefix}OvlErrorRate <float=unset>
  Overlaps above this error rate are not computed.
  * `corOvlErrorRate` applies to overlaps generated for correcting reads;
  * `obtOvlErrorRate` applied to overlaps generated for trimming reads;
  * `utgOvlErrorRate` applies to overlaps generated for assembling reads.
  These limits apply to the 'ovl' overlap algorithm and when alignments are computed for mhap overlaps with :ref:`mhapReAlign <mhapReAlign>`.

{prefix}OvlFrequentMers <string=undefined>
  Do not seed overlaps with these kmers (fasta format).

{prefix}OvlHashBits <integer=unset>
  Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per corOvlHashBlockLength.

{prefix}OvlHashBlockLength <integer=unset>
  Amount of sequence (bp to load into the overlap hash table.

{prefix}OvlHashLoad <integer=unset>
  Maximum hash table load.  If set too high, table lookups are inefficent; if too low, search overhead dominates run time.

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
  Number of reads per 1GB block.  Memory * size is loaded into memory per job.

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

Overlapper Configuration, mhap Algorithm
----------------------------------------

{prefix}MMapBlockSize <integer=unset>
  Number of reads per 1GB block.  Memory * size is loaded into memory per job.

{prefix}MMapMerSize <integer=unset>
  K-mer size for seeds in minimap.

Overlap Store
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The overlap algorithms return overlaps in an arbitrary order.  The correction, trimming and assembly
algorithms usually need to know all overlaps for a single read.  The overlap store duplicates each
overlap, sorts them by the first ID, and stores them for quick retrieval of all overlaps for a
single read.

ovsMemory <float>
  How much memory, in gigabytes, to use for constructing overlap stores.  Must be at least 256m or 0.25g.

ovsMethod <string="sequential">
  Two construction algorithms are supported.  One uses a single data stream, and is faster for small
  and moderate size assemblies.  The other uses parallel data streams and can be faster (depending
  on your network disk bandwitdh) for moderate and large assemblies.

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



.. _grid-engine:

Grid Engine Support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Canu directly supports most common grid scheduling systems.  Under normal use, Canu will query the
system for grid support, congigure itself for the machines available in the grid, then submit itself
to the grid for execution.  The Canu pipeline is a series of about a dozen steps that alternate
between embarassingly parallel computations (e.g., overlap computation) and sequential bookkeeping
steps (e.g., checking if all overlap jobs finished).  This is entirely managed by Canu.

Canu has first class support for the various schedulers derived from Sun Grid Engine (Univa, Son of
Grid Engine) and the Simple Linux Utility for Resource Management (SLURM), meaning that the
devlopers have direct access to these systems.  Platform Computing's Load Sharing Facility (LSF) and
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

For SGE grids, two options are sometimes necessary to tell canu about pecularities of your grid:
``gridEngineThreadsOption`` describes how to request multiple cores, and ``gridEngineMemoryOption``
describes how to request memory.  Usually, canu can figure out how to do this, but sometimes it
reports an error such as::

 -- WARNING:  Couldn't determine the SGE parallel environment to run multi-threaded codes.
 --           Valid choices are (pick one and supply it to canu):
 --             gridEngineThreadsOption="-pe make THREADS"
 --             gridEngineThreadsOption="-pe make-dedicated THREADS"
 --             gridEngineThreadsOption="-pe mpich-rr THREADS"
 --             gridEngineThreadsOption="-pe openmpi-fill THREADS"
 --             gridEngineThreadsOption="-pe smp THREADS"
 --             gridEngineThreadsOption="-pe thread THREADS"

or::

 -- WARNING:  Couldn't determine the SGE resource to request memory.
 --           Valid choices are (pick one and supply it to canu):
 --             gridEngineMemoryOption="-l h_vmem=MEMORY"
 --             gridEngineMemoryOption="-l mem_free=MEMORY"

If you get such a message, just add the appropriate line to your canu command line.  Both options
will replace the uppercase text (THREADS or MEMORY) with the value canu wants when the job is
submitted.  For ``gridEngineMemoryOption``, any number of ``-l`` options can be supplied; we could
use ``gridEngineMemoryOption="-l h_vmem=MEMORY -l mem_free=MEMORY"`` to request both ``h_vmem`` and
``mem_free`` memory.

.. _grid-options:

Grid Options
------------

To run on the grid, each stage needs to be configured - to tell the grid how many cores it will use and how much memory it needs.
Some support for this is automagic (for example, overlapInCore and mhap know how to do this), others need to be manually configured.
Yes, it's a problem, and yes, we want to fix it.

The gridOptions* parameters supply grid-specific opitons to the grid submission command.

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
assebmled, the type of processing desired, or the amount of comput resources available.  Overlap

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

+-------+-----------+----------------------------------------+
|    Prefix         | Algorithm                              |
+=======+===========+========================================+
| cor   |           | | Overlap generation using the         |
+-------|           | | 'mhap' algorithm for                 |
| obt   | mhap      | | 'cor'=correction,, 'obt'=trimming    |
+-------|           | | or 'utg'=assembly.                   |
| utg   |           |                                        |
+-------+-----------+----------------------------------------+
| cor   |           | | Overlap generation using the         |
+-------|           | | 'minimap' algorithm for              |
| obt   | mmap      | | 'cor'=correction,, 'obt'=trimming    |
+-------|           | | or 'utg'=assembly.                   |
| utg   |           |                                        |
+-------+-----------+----------------------------------------+
| cor   |           | | Overlap generation using the         |
+-------|           | | 'overlapInCore' algorithm for        |
| obt   | ovl       | | 'cor'=correction,, 'obt'=trimming    |
+-------|           | | or 'utg'=assembly.                   |
| utg   |           |                                        |
+-------+-----------+----------------------------------------+
|       | ovb       | Parallel overlap store bucketizing     |
+-------+-----------+----------------------------------------+
|       | ovs       | Parallel overlap store bucket sorting  |
+-------+-----------+----------------------------------------+
|       | cor       | Read correction                        |
+-------+-----------+----------------------------------------+
|       | red       | Error detection in reads               |
+-------+-----------+----------------------------------------+
|       | oea       | Error adjustment in overlaps           |
+-------+-----------+----------------------------------------+
|       | bat       | Unitig/contig construction             |
+-------+-----------+----------------------------------------+
|       | cns       | Unitig/contig consensus                |
+-------+-----------+----------------------------------------+

For example, 'mhapMemory` would set the memory limit for computing overlaps with the mhap algorithm;
'cormhapMemory' would set the memory limit only when mhap is used for generating overlaps used for
correction.

The 'minMemory', 'maxMemory', 'minThreads' and 'maxThreads' options will apply to all jobs, and
can be used to artifically limit canu to a portion of the current machine.  In the overlapper
example above, setting maxThreads=4 would result in two concurrent jobs instead of four.


Overlap Error Adjustment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

red = Read Error Detection
oea = Overlap Error Adjustment

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
  Don't make a paritition with fewer than N reads

cnsMaxCoverage
  Limit unitig consensus to at most this coverage.
 
.. _cnsErrorRate:

cnsErrorRate
  Inform the consensus genration algorithm of the amount of difference it should expect in a
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

corErrorRate <integer=unset>
  Do not use overlaps with error rate higher than this (estimated error rate for `mhap` and `minimap` overlaps).

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

By default, canu will split the final output into three files:

asm.contigs.fasta
   Everything which could be assembled and is part of the primary assembly, including both unique and repetitive elements.  Each contig has several flags included on the fasta def line::

asm.bubbles.fasta
   alternate paths in the graph which could not be merged into the primary assembly.

asm.unassembled.fasta
   reads/tigs which could not be incorporated into the primary or bubble assemblies.

It is possible for tigs comprised of multiple reads to end up in asm.unassembled.fasta. The default filtering eliminates anything with < 2 reads, shorter than 1000bp, or comprised of mostly a single sequence (>75%). The filtering is controlled by the contigFilter parameter which takes 5 values.

::

   contigFilter
     minReads
     minLength
     singleReadSpan
     lowCovSpan
     lowCovDepth

The default filtering is "2 1000 0.75 0.75 2". If you are assembling amplified data or viral data, it is possible your assembly will be flagged as unassembled. In those cases, you can turn off the filtering with the parameters

::

   contigFilter="2 1000 1.0 1.0 2"
