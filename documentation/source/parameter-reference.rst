
.. _parameter-reference:

Canu Parameter Reference
========================

To get the most up-to-date options, run

   canu -options

The default values below will vary based on the input data type and genome size.

.. _genomeSize:

Global Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The catch all category.

errorRate <float=0.01>
  The expected error rate, as fraction error, in the corrected reads, set by default based on data type, typically not changed by the user.

genomeSize <float=unset>
  An estimate of the size of the genome.  Common suffices are allowed, for example, 3.7m or 2.8g. Required.

canuIteration <internal parameter, do not use>
  Which parallel iteration is being attempted.
canuIterationMax <integer=2>
  How many parallel iterations to try.  Ideally, the parallel jobs, run under grid control, would all finish successfully on the first try.
  Sometimes, jobs fail due to other jobs exhausting resources (memory), or by the node itself failing.  In this case, canu will launch the jobs
  again.  This parameter controls how many times it tries.


Process Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

showNext <boolean=false>
  Report the first major command that would be run, but don't run it.  Processing to get to that
  command, for example, checking the output of the previous command or preparing inputs for the
  next command, is still performed.

stopBefore <string=undefined>
  Stop processing just before this stage would execute.  The stage is configured, and the
  command logged to the standard output, before stopping.  For grid-based stages, e.g., overlapper,
  the grid submit command is reported.

  If the stage has finished successfully, it will not stop.

  Only one stage may be supplied to stopBefore.

  Valid stages to stopBefore are:
    - gatekeeper
    - meryl
    - overlap
    - correction
    - overlapErrorAdjustment
    - trimReads
    - splitReads
    - unitig
    - consensusConfigure
    - consensus
    - output

  The default value is 'undef'.

stopAfter <string=undefined>
  Stop processing after this stage completes.

  The stage will typically stop BEFORE a summary of its processing is reported in the canu chatter.

  Valid stages top stopAfter are:
    - gatekeeper
    - meryl
    - overlapConfigure
    - overlap
    - correction
    - unitig
    - consensusConfigure
    - consensus
    - consensusLoad
    - consensusFilter


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

{prefix}OvlErrorRate <float=unset>
  Overlaps above this error rate are not computed.

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

{prefix}ReAlign <boolean=false>
  Compute actual alignments from mhap overlaps; 'raw' from mhap output;
  uses either obtErrorRate or ovlErrorRate, depending on which overlaps are computed)

{prefix}MhapSensitivity <string="normal">
  Coarse sensitivity level: 'normal' or 'high' or 'fast'.

Overlapper Configuration, mhap Algorithm
----------------------------------------

{prefix}MMapBlockSize <integer=unset>
  Number of reads per 1GB block.  Memory * size is loaded into memory per job.

{prefix}MMapMerSize <integer=unset>
  K-mer size for seeds in minimap.

{prefix}ReAlign <boolean=false>
  Compute actual alignments from minimap overlaps; 'raw' from mhap output;
  uses either obtErrorRate or ovlErrorRate, depending on which overlaps are computed)

Overlap Store
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The overlap algorithms return overlaps in an arbitrary order.  The correction, trimming and assembly
algorithms usually need to know all overlaps for a single read.  The overlap store duplicates each
overlap, sorts them by the first ID, and stores them for quick retrieval of all overlaps for a
single read.

ovsMemory <integer=2>
  How much memory, in gigabytes, to use for constructing overlap stores.

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

obtErrorRate <float=unset>
  Stringency of overlaps to use for trimming

trimReadsOverlap <integer=1>
  Minimum overlap between evidence to make contiguous trim.

trimReadsCoverage <integer=1>
  Minimum depth of evidence to retain bases.



.. _grid-engine:

Grid Engine Support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Canu supports Sun/Open/Univa Grid Engine (SGE), Portable Batch System (PBS), Platform Computing's Load
Sharing Facility (LSF), and the Simple Linux Utility for Resource Management (SLURM).  Most of the compute
intensive stages can run under grid control.

The useGrid* options control which algorithms run in parallel on the grid.

useGrid <boolean=true>
  Master control.  If 'false', no algorithms will run under grid control.  Does not change the value of the other useGrid options.

  If 'remote', jobs are configured for grid execution, but not submitted.  A message, with commands to launch the job, is reported and canu halts execution.

  Note that the host used to run canu for 'remote' execution must know about the grid, that is, it must be able to submit jobs to the grid.

There are many options for configuring a new grid ('gridEngine*') and for configuring how canu
configures its computes to run under grid control ('gridOptions*').  The grid engine to use is
specified with the 'gridEngine' option.

gridEngine <string>
  Which grid engine to use.  Auto-detected.  Possible choices are 'sge', 'pbs', 'lsf' or 'slurm'.

  NOTE: 'lsf' support is untested.

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

WHERE IS OBT??


Algorithm Execution Method
--------------------------

Each of the high compute stages can be computed either on a grid or in parallel on the local machine.
Most algorithms will respect a given maximum memory usage.
Most algorithms can support more than a single thread of computation.
When the grid engine is not used, more than one task can be run at a time.

BUG:  Local execution doesn't pay attention to memory option.

For execution locally, three parameters describe the task:

{prefix}Concurrency <integer=unset>
  Set the number of tasks that can run at the same time, when running without grid support.

  Available prefixes are:
    - master
    - cns
    - cor
    - cormhap
    - obtmhap
    - utgmhap
    - corovl
    - obtovl
    - utgovl
    - cormmap
    - obtmmap
    - utgmmap
    - oea
    - ovb
    - ovs
    - red

{prefix}Threads <integer=unset>
  Set the number of compute threads used per task.

  Available prefixes are:
    - master
    - bat
    - cns
    - cor
    - cormhap
    - obtmhap
    - utgmhap
    - corovl
    - obtovl
    - utgovl
    - cormmap
    - obtmmap
    - utgmmap
    - ovb
    - ovs
    - red
    - oea

{prefix}Memory <integer=unset>
  Set the amount of memory, in GB, to use for each job in a task.

  Available prefixes are:
    - master
    - bat
    - ovb
    - ovs
    - cns
    - cor
    - cormhap
    - obtmhap
    - utgmhap
    - corovl
    - obtovl
    - utgovl
    - cormmap
    - obtmmap
    - utgmmap
    - red
    - oea

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
 
cnsErrorRate
  Possibly unused.


Read Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in Canu is to find high-error overlaps and generate corrected sequences for subsequent assembly. This is currently the fastest step in Canu. By default, only the longest 40X of data (based on the specified genome size) is used for correction. Typically, some reads are trimmed during correction due to being chimeric or having erroneous sequence, resulting in a loss of 20-25% (30X output). You can force correction to be non-lossy by setting 

::

   corMinCoverage=0

In which case the corrected reads output will be the same length as the input data, keeping any high-error unsupported bases. Canu will trim these in downstream steps before assembly.

If you have a dataset with uneven coverage or small plasmids, correcting the longest 40X may not give you sufficient coverage of your genome/plasmid. In these cases, you can set 

::

   corOutCoverage=400

Or any large value greater than your total input coverage which will correct and assemble all input data, at the expense of runtime.

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

corOutCoverage <integer=40>
  Only correct the longest reads up to this coverage; default 40

corFilter <string="expensive">
  Method to filter short reads from correction; 'quick' or 'expensive' or 'none'

falconSense
  Path to fc_consensus.py or falcon_sense.bin

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
