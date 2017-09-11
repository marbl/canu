
.. _faq:

Canu FAQ
========


.. contents::
  :local:


What resources does Canu require for a bacterial genome assembly? A mammalian assembly?
-------------------------------------
    Canu will detect available resources and configure itself to run efficiently using those
    resources.  It will request resources, for example, the number of compute threads to use, Based
    on the genome size being assembled. It will fail to even start if it feels there are
    insufficient resources available.
    
    A typical bacterial genome can be assembled with 8GB memory in a few CPU hours - around an hour
    on 8 cores.  It is possible, but not allowed by default, to run with only 4GB memory.

    A well-behaved large genome, such as human or other mammals, can be assembled in 10,000 to
    25,000 CPU hours, depending on coverage.  A grid environment is strongly recommended, with at
    least 16GB available on each compute node, and one node with at least 64GB memory.  You should
    plan on having 3TB free disk space, much more for highly repetitive genomes.

    Our compute nodes have 48 compute threads and 128GB memory, with a few larger nodes with up to
    1TB memory.  We develop and test (mostly bacteria, yeast and drosophila) on laptops and desktops
    with 4 to 12 compute threads and 16GB to 64GB memory.

    
How do I run Canu on my SLURM / SGE / PBS / LSF / Torque system?
-------------------------------------
    Canu will detect and configure itself to use on most grids. You can supply your own grid
    options, such as a partition on SLURM or an account code on SGE, with ``gridOptions="<your
    options list>"`` which will passed to every job submitted by Canu.  Similar options exist for
    every stage of Canu, which could be used to, for example, restrict overlapping to a specific
    partition or queue.

    To disable grid support and run only on the local machine, specify ``useGrid=false``
    
My run stopped with the error ``'Failed to submit batch jobs'``
-------------------------------------

    The grid you run on must allow compute nodes to submit jobs. This means that if you are on a compute host, ``qsub/bsub/sbatch/etc`` must be available and working. You can test this by starting an interactive compute session and running the submit command manually (e.g. ``qsub`` on SGE, ``bsub`` on LSF, ``sbatch`` on SLURM). 
    
    If this is not the case, Canu **WILL NOT** work on your grid. You must then set ``useGrid=false`` and run on a single machine. Alternatively, you can run Canu with ``useGrid=remote`` which will stop at every submit command, list what should be submitted. You then submit these jobs manually, wait for them to complete, and run the Canu command again. This is a manual process but currently the only workaround for grids without submit support on the compute nodes.


What parameters should I use for my reads?
-------------------------------------
    Canu is designed to be universal on a large range of PacBio (C2, P4-C2, P5-C3, P6-C4) and Oxford Nanopore
    (R6 through R9) data.  Assembly quality and/or efficiency can be enhanced for specific datatypes:
    
    **Nanopore R7 1D** and **Low Identity Reads**
       With R7 1D sequencing data, and generally for any raw reads lower than 80% identity, five to
       ten rounds of error correction are helpful::

         canu -p r1 -d r1 -correct corOutCoverage=500 corMinCoverage=0 corMhapSensitivity=high -nanopore-raw your_reads.fasta
         canu -p r2 -d r2 -correct corOutCoverage=500 corMinCoverage=0 corMhapSensitivity=high -nanopore-raw r1/r1.correctedReads.fasta.gz
         canu -p r3 -d r3 -correct corOutCoverage=500 corMinCoverage=0 corMhapSensitivity=high -nanopore-raw r2/r2.correctedReads.fasta.gz
         canu -p r4 -d r4 -correct corOutCoverage=500 corMinCoverage=0 corMhapSensitivity=high -nanopore-raw r3/r3.correctedReads.fasta.gz
         canu -p r5 -d r5 -correct corOutCoverage=500 corMinCoverage=0 corMhapSensitivity=high -nanopore-raw r4/r4.correctedReads.fasta.gz

       Then assemble the output of the last round, allowing up to 30% difference in overlaps::

         canu -p asm -d asm correctedErrorRate=0.3 utgGraphDeviation=50 -nanopore-corrected r5/r5.correctedReads.fasta.gz

    **Nanopore R7 2D** and **Nanopore R9 1D**
      The defaults were designed with these datasets in mind so they should work.

    **Nanopore R9 2D** and **PacBio P6**
       Slightly decrease the maximum allowed difference in overlaps from the default of 14.4% to 12.0%
       with ``correctedErrorRate=0.120``

    **Early PacBio Sequel**
       Based on exactly one publically released *A. thaliana* `dataset
       <http://www.pacb.com/blog/sequel-system-data-release-arabidopsis-dataset-genome-assembly/>`_,
       slightly decrease the maximum allowed difference from the default of 4.5% to 4.0% with
       ``correctedErrorRate=0.040 corMhapSensitivity=normal``.  For recent Sequel data, the defaults
       seem to be appropriate.
       
   **Nanopore R9 large genomes**
       Due to some systematic errors, the identity estimate used by Canu for correction can be an
       over-estimate of true error, inflating runtime. For recent large genomes (>1gbp) with more
       than 30x coverage, we've used ``'corMhapOptions=--threshold 0.8 --num-hashes
       512 --ordered-sketch-size 1000 --ordered-kmer-size 14'``. This is not needed for below 30x
       coverage.


My assembly continuity is not good, how can I improve it?
-------------------------------------
    The most important determinant for assembly quality is sequence length, followed by the repeat
    complexity/heterozygosity of your sample.  The first thing to check is the amount of corrected
    bases output by the correction step.  This is logged in the stdout of Canu or in
    canu-scripts/canu.*.out if you are running in a grid environment. For example on `a
    haploid H. sapiens <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN02744161>`_ sample:
    
    ::
    
       -- BEGIN TRIMMING
       --
       ...
       -- In gatekeeper store 'chm1/trimming/asm.gkpStore':
       --   Found 5459105 reads.
       --   Found 91697412754 bases (29.57 times coverage).
       ...

   Canu tries to correct the longest 40X of data. Some loss is normal but having output coverage
   below 20-25X is a sign that correction did not work well (assuming you have more input coverage
   than that). If that is the case, re-running with ``corMhapSensitivity=normal`` if you have >50X
   or ``corMhapSensitivity=high corMinCoverage=0`` otherwise can help. You can also increase the
   target coverage to correct ``corOutCoverage=100`` to get more correct sequences for assembly. If
   there are sufficient corrected reads, the poor assembly is likely due to either repeats in the
   genome being greater than read lengths or a high heterozygosity in the sample. Stay tuned for mor
   information on tuning unitigging in those instances.


.. _tweak:

What parameters can I tweak?
-------------------------------------
    For all stages:

    - ``rawErrorRate`` is the maximum expected difference in an alignment of two _uncorrected_
      reads.  It is a meta-parameter that sets other parameters.

    - ``correctedErrorRate`` is the maximum expected difference in an alignment of two _corrected_
      reads.  It is a meta-parameter that sets other parameters.  (If you're used to the
      ``errorRate`` parameter, multiply that by 3 and use it here.)

    - ``minReadLength`` and ``minOverlapLength``.  The defaults are to discard reads shorter than
      1000bp and to not look for overlaps shorter than 500bp.  Increasing ``minReadLength`` can
      improve run time, and increasing ``minOverlapLength`` can improve assembly quality by removing
      false overlaps.  However, increasing either too much will quickly degrade assemblies by either
      omitting valuable reads or missing true overlaps.

    For correction:

    - ``corOutCoverage`` controls how much coverage in corrected reads is generated.  The default is
      to target 40X, but, for various reasons, this results in 30X to 35X of reads being generated.

    - ``corMinCoverage``, loosely, controls the quality of the corrected reads.  It is the coverage
      in evidence reads that is needed before a (portion of a) corrected read is reported.
      Corrected reads are generated as a consensus of other reads; this is just the minimum ocverage
      needed for the consensus sequence to be reported.  The default is based on input read
      coverage: 0x coverage for less than 30X input coverage, and 4x coverage for more than that.

    For assembly:

    - ``utgOvlErrorRate`` is essientially a speed optimization.  Overlaps above this error rate are
      not computed.  Setting it too high generally just wastes compute time, while setting it too
      low will degrade assemblies by missing true overlaps between lower quality reads.

    - ``utgGraphDeviation`` and ``utgRepeatDeviation`` what quality of overlaps are used in contig
      construction or in breaking contigs at false repeat joins, respectively.  Both are in terms of
      a deviation from the mean error rate in the longest overlaps.

    - ``utgRepeatConfusedBP`` controls how similar a true overlap (between two reads in the same
      contig) and a false overlap (between two reads in different contigs) need to be before the
      contig is split.  When this occurs, it isn't clear which overlap is 'true' - the longer one or
      the slightly shorter one - and the contig is split to avoid misassemblies.

    For polyploid genomes:

        Generally, there's a couple of ways of dealing with the ploidy. 
    
        1) **Avoid collapsing the genome** so you end up with double (assuming diploid) the genome
           size as long as your divergence is above about 2% (for PacBio data). Below this
           divergence, you'd end up collapsing the variations. We've used the following parameters
           for polyploid populations (PacBio data):

           ``corOutCoverage=200 correctedErrorRate=0.040 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50"``
    
           This will output more corrected reads (than the default 40x). The latter option will be
           more conservative at picking the error rate to use for the assembly to try to maintain
           haplotype separation. If it works, you'll end up with an assembly >= 2x your haploid
           genome size. Post-processing using gene information or other synteny information is
           required to remove redunancy from this assembly.

        2) **Smash haplotypes together** and then do phasing using another approach (like HapCUT2 or
           whatshap or others). In that case you want to do the opposite, increase the error rates
           used for finding overlaps:
   
           ``corOutCoverage=200 ovlErrorRate=0.15 obtErrorRate=0.15``

           Error rates for trimming (``obtErrorRate``) and assembling (``batErrorRate``) can usually
           be left as is.  When trimming, reads will be trimmed using other reads in the same
           chromosome (and probably some reads from other chromosomes).  When assembling, overlaps
           well outside the observed error rate distribution are discarded.

    For low coverage:

     - For less than 30X coverage, increase the alllowed difference in overlaps from 4.5% to 7.5%
       (or more) with ``correctedErrorRate=0.075``, to adjust for inferior read correction.  Canu
       will automatically reduce ``corMinCoverage`` to zero to correct as many reads as possible.

    For high coverage:

     - For more than 60X coverage, decrease the allowed difference in overlaps from 4.5% to 4.0%
       with ``correctedErrorRate=0.040``, so that only the better corrected reads are used.  This is
       primarily an optimization for speed and generally does not change assembly continuity.


My asm.contigs.fasta is empty, why?
-------------------------------------
    Canu creates three assembled sequence :ref:`output files <outputs>`: ``<prefix>.contigs.fasta``,
    ``<prefix>.unitigs.fasta``, and ``<prefix>.unassembled.fasta``, where contigs are the primary
    output, unitigs are the primary output split at alternate paths,
    and unassembled are the leftover pieces.

    The :ref:`contigFilter <contigFilter>` parameter sets several parameters that control how small
    or low coverage initial contigs are handled.  By default, initial contigs with more than 50% of
    the length at less than 5X coverage will be classified as 'unassembled' and removed from the
    assembly, that is, ``contigFilter="2 0 1.0 0.5 5"``.  The filtering can be disabled by changing
    the last number from '5' to '0' (meaning, filter if 50% of the contig is less than 0X coverage).


Why is my assembly is missing my favorite short plasmid?
-------------------------------------
    Only the longest 40X of data (based on the specified genome size) is used for
    correction.  Datasets with uneven coverage or small plasmids can fail to generate enough
    corrected reads to give enough coverage for assembly, resulting in gaps in the genome or even no
    reads for small plasmids.  Set ``corOutCoverage=1000`` (or any value greater than your total input
    coverage) to correct all input data.

    An alternate approach is to correct all reads (``-correct corOutCoverage=1000``) then assemble
    40X of reads picked at random from the ``<prefix>.correctedReads.fasta.gz`` output.


Why do I get less corrected read data than I asked for?
-------------------------------------
    Some reads are trimmed during correction due to being chimeric or because there wasn't enough
    evidence to generate a quality corrected sequence.  Typically, this results in a 25% loss.
    Setting ``corMinCoverage=0`` will report all bases, even low those of low quality.  Canu will
    trim these in its 'trimming' phase before assembly.


What is the minimum coverage required to run Canu?
-------------------------------------
    For eukaryotic genomes, coverage more than 20X is enough to outperform current hybrid methods.


My genome is AT (or GC) rich, do I need to adjust parameters?  What about highly repetitive genomes?
-------------------------------------
   On bacterial genomes, no adjustment of parameters is (usually) needed.  See the next question.

   On repetitive genomes with with a significantly skewed AT/GC ratio, the Jaccard estimate used by
   MHAP is biased.  Setting ``corMaxEvidenceErate=0.15`` is sufficient to correct for the bias in
   our testing.

   In general, with high coverage repetitive genomes (such as plants) it can be beneficial to set
   the above parameter anyway, as it will eliminate repetitive matches, speed up the assembly, and
   sometime improve unitigs.


How can I send data to you?
-------------------------------------
   FTP to ftp://ftp.cbcb.umd.edu/incoming/sergek.  This is a write-only location that only the Canu
   developers can see.
   
   Here is a quick walk-through using a command-line ftp client (should be available on most Linux and OSX installations). Say we want to transfer a file named ``reads.fastq``. First, run ``ftp ftp.cbcb.umd.edu``, specify ``anonymous`` as the user name and hit return for password (blank). Then:
   
   .. code-block:: 
   
      cd incoming/sergek
      put reads.fastq
      quit

   That's it, you won't be able to see the file but we can download it.
