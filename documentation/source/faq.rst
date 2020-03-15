
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
    25,000 CPU hours, depending on coverage from raw data.  A well-behaved large genome, such as human or other mammals with HiFi data can be assembled in 1,000 to 2,000 CPU hours. A grid environment is strongly recommended, with at
    least 16GB available on each compute node, and one node with at least 64GB memory.  You should
    plan on having 3TB free disk space (200 GB for HiFi), much more for highly repetitive genomes.

    Our compute nodes have 48 compute threads and 128GB memory, with a few larger nodes with up to
    1TB memory.  We develop and test (mostly bacteria, yeast and drosophila) on laptops and desktops
    with 4 to 12 compute threads and 16GB to 64GB memory.


How do I run Canu on my SLURM / SGE / PBS / LSF / Torque system?
-------------------------------------
    Canu will detect and configure itself to use on most grids. Canu will NOT request explicit time limits or
    queues/partitions. You can supply your own grid options, such as a partition on SLURM, an account code 
    on SGE, and/or time limits with ``gridOptions="<your options list>"`` which will passed to every job 
    submitted by Canu.  Similar options exist for every stage of Canu, which could be used to, for example, 
    restrict overlapping to a specific partition or queue.

    To disable grid support and run only on the local machine, specify ``useGrid=false``

    It is possible to limit the number of grid jobs running at the same time, but this isn't
    directly supported by Canu.  The various :ref:`gridOptions <grid-options>` parameters
    can pass grid-specific parameters to the submit commands used; see
    `Issue #756 <https://github.com/marbl/canu/issues/756>`_ for Slurm and SGE examples.


My run stopped with the error ``'Mhap precompute jobs failed'``
-------------------------------------

    Several package managers make a mess of the installation causing this error (conda and ubuntu in particular). Package managers don't add much benefit to a tool like Canu which is distributed as pre-compiled binaries compatible with most systems so our recommended installation method is downloading a binary release. Try running the assembly from scratch using our release distribution and if you continue to encounter errors, submit an issue.

My run stopped with the error ``'Failed to submit batch jobs'``
-------------------------------------

    The grid you run on must allow compute nodes to submit jobs. This means that if you are on a
    compute host, ``qsub/bsub/sbatch/etc`` must be available and working. You can test this by
    starting an interactive compute session and running the submit command manually (e.g. ``qsub``
    on SGE, ``bsub`` on LSF, ``sbatch`` on SLURM).

    If this is not the case, Canu **WILL NOT** work on your grid. You must then set
    ``useGrid=false`` and run on a single machine. Alternatively, you can run Canu with
    ``useGrid=remote`` which will stop at every submit command, list what should be submitted. You
    then submit these jobs manually, wait for them to complete, and run the Canu command again. This
    is a manual process but currently the only workaround for grids without submit support on the
    compute nodes.


My run of Canu was killed by the sysadmin; the power going out; my cat stepping on the power button; et cetera.  Is it safe to restart?  How do I restart?
-------------------------------------

    Yes, perfectly safe!  It's actually how Canu runs normally: each time Canu starts, it examines
    the state of the assembly to decide what it should do next.  For example, if six overlap tasks
    have no results, it'll run just those six tasks.

    This also means that if you want to redo some step, just remove those results from the assembly
    directory.  Some care needs to be taken to make sure results computed after those are also
    removed.

    Short answer: just rerun the _exact_ same command as before.  It'll do the right thing.


My genome size and assembly size are different, help!
-------------------------------------
    The difference could be due to a heterozygous genome where the assembly separated some loci. It could also be because the previous estimate is incorrect. We typically use two analyses to see what happened. First, a `BUSCO <https://busco.ezlab.org>`_ analysis will indicate duplicated genes. For example this assembly::

      INFO	C:98.5%[S:97.9%,D:0.6%],F:1.0%,M:0.5%,n:2799
      INFO	2756 Complete BUSCOs (C)
      INFO	2740 Complete and single-copy BUSCOs (S)
      INFO	16 Complete and duplicated BUSCOs (D)
    
    does not have much duplication but this assembly::
    
      INFO	C:97.6%[S:15.8%,D:81.8%],F:0.9%,M:1.5%,n:2799
      INFO	2732 Complete BUSCOs (C)
      INFO	443 Complete and single-copy BUSCOs (S)
      INFO	2289 Complete and duplicated BUSCOs (D)
    
    does. We have had success using `purge_dups <https://github.com/dfguan/purge_dups>`_ to remove duplication. Purge dups will also generate a coverage histogram which will usually have two peaks when assemblies have separated some loci, make sure to inspect it to make sure the cutoffs selected are reasonable. 

What parameters should I use for my reads?
-------------------------------------
    Canu is designed to be universal on a large range of PacBio CLR, PacBio HiFi, Oxford
    Nanopore (R6 through R10) data.  Assembly quality and/or efficiency can be enhanced for specific
    datatypes:

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
      The defaults were designed with these datasets in mind so they should work. Having very high
      coverage or very long Nanopore reads can slow down the assembly significantly. You can try the
      ``-fast`` option which is much faster but may produce less
      contiguous assemblies on large genomes.

    **Nanopore flip-flop R9.4 or R10.3**
       Based on a human dataset, the flip-flop basecaller reduces both the raw read error rate and the residual error rate remaining after Canu read correction. For this reason you can reduce the error tolerated by Canu. If you have over 30x coverage add the options: ``'corMhapOptions=--threshold 0.8 --ordered-sketch-size 1000 --ordered-kmer-size 14' correctedErrorRate=0.105``. This is primarily a speed optimization so you can use defaults, especially if your genome's accuracy is not improved by the flip-flop caller.

    **PacBio Sequel**
       Based on an *A. thaliana* `dataset
       <http://www.pacb.com/blog/sequel-system-data-release-arabidopsis-dataset-genome-assembly/>`_,
       and a few more recent mammalian genomes, slightly increase the maximum allowed difference from the default of 4.5% to 8.5% with
       ``correctedErrorRate=0.085 corMhapSensitivity=normal``.
      Only add the second parameter (``corMhapSensivity=normal``) if you have >50x coverage.

    **PacBio Sequel II**
       The defaults for PacBio should work on this data. You could speed up the assembly by decreasing the error rate from the default, especially if you have high (>50x) coverage via ``correctedErrorRate=0.035 utgOvlErrorRate=0.065 trimReadsCoverage=2 trimReadsOverlap=500``

    **PacBio HiFi**
       The defaults for -pacbio-hifi should work on this data. There is still some variation in data quality between samples. If you have poor continuity, it may be because the data is lower quality than expected. Try running the assembly with ``-trim-assemble``` or with ``batOptions="-eg 0.01 -sb 0.01 -dg 6 -db 6 -dr 1 -ca 50 -cp 5”``. You will likely get a genome size larger than you expect, due to separation of alleles. See `My genome size and assembly size are different, help!`_ for details on how to remove this duplication.


Can I assemble RNA sequence data?
-------------------------------------
    Canu will likely mis-assemble, or completely fail to assemble, RNA data.  It will do a
    reasonable job at generating corrected reads though.  Reads are corrected using (local) best
    alignments to other reads, and alignments between different isoforms are usually obviously not
    'best'.  Just like with DNA sequences, similar isoforms can get 'mixed' together.  We've heard
    of reasonable success from users, but do not have any parameter suggestions to make.

    Note that Canu will silently translate 'U' bases to 'T' bases on input, but **NOT** translate
    the output bases back to 'U'.
    
Can I assemble amplicon sequence data?
-------------------------------------   
    In short, yes. Typically these have very high coverage so we recommend randomly downsampling (``'readSamplingCoverage=100'``) and turning off filtering of short contigs ``contigFilter="2 0 1.0 0.5 0"``.

My assembly is running out of space, is too slow?
-------------------------------------
    We don't have a good way to estimate of disk space used for the assembly. It varies with genome size, repeat content, and sequencing depth. A human genome sequenced with PacBio or Nanopore at 40-50x typically requires 1-2TB of space at the peak. Plants, unfortunately, seem to want a lot of space. 10TB is a reasonable guess. We've seen it as bad as 20TB on some very repetitive genomes.
    
    The most common cause of high disk usage is a very repetitive or large genome. There are some parameters you can tweak to both reduce disk space and speed up the run. Try adding the options ``corMhapFilterThreshold=0.0000000002 corMhapOptions="--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50" mhapMemory=60g mhapBlockSize=500 ovlMerDistinct=0.975``. This will suppress repeats more than the default settings and speed up both correction and assembly.
    
    It is also possible to clean up some intermediate outputs before the assembly is complete to save space. If you already have a ```*.ovlStore.BUILDING/1-bucketize.successs`` file in your current step (e.g. ``correct```), you can clean up the files under ``1-overlapper/blocks``. You can also remove the ovlStore for the previous step if you have its output (e.g. if you have ``asm.trimmedReads.fasta.gz``, you can remove ``trimming/asm.ovlStore``). 

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
      Corrected reads are generated as a consensus of other reads; this is just the minimum coverage
      needed for the consensus sequence to be reported.  The default is based on input read
      coverage: 0x coverage for less than 30X input coverage, and 4x coverage for more than that.

    For assembly:

    - ``utgOvlErrorRate`` is essentially a speed optimization.  Overlaps above this error rate are
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

           ``corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50"``

           This will output more corrected reads (than the default 40x). The latter option will be
           more conservative at picking the error rate to use for the assembly to try to maintain
           haplotype separation. If it works, you'll end up with an assembly >= 2x your haploid
           genome size. Post-processing using gene information or other synteny information is
           required to remove redundancy from this assembly.

        2) **Smash haplotypes together** and then do phasing using another approach (like HapCUT2 or
           whatshap or others). In that case you want to do the opposite, increase the error rates
           used for finding overlaps:

           ``corOutCoverage=200 correctedErrorRate=0.15``

           When trimming, reads will be trimmed using other reads in the same
           chromosome (and probably some reads from other chromosomes).  When assembling, overlaps
           well outside the observed error rate distribution are discarded.
           
         We strongly recommend option 1 which will lead to a larger than expected genome size. See `My genome size and assembly size are different, help!`_ for details on how to remove this duplication.

    For metagenomes:

        The basic idea is to use all data for assembly rather than just the longest as default. The
        parameters we've used recently are:

          ``corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=200``

    For low coverage:

     - For less than 30X coverage, increase the alllowed difference in overlaps by a few percent
       (from 4.5% to 8.5% (or more) with ``correctedErrorRate=0.105`` for PacBio and from 14.4% to
       16% (or more) with ``correctedErrorRate=0.16`` for Nanopore), to adjust for inferior read
       correction.  Canu will automatically reduce ``corMinCoverage`` to zero to correct as many
       reads as possible.

    For high coverage:

     - For more than 60X coverage, decrease the allowed difference in overlaps (from 4.5% to 4.0%
       with ``correctedErrorRate=0.040`` for PacBio, from 14.4% to 12% with
       ``correctedErrorRate=0.12`` for Nanopore), so that only the better corrected reads are used.
       This is primarily an optimization for speed and generally does not change assembly
       continuity.


My asm.contigs.fasta is empty, why?
-------------------------------------
    Canu creates three assembled sequence :ref:`output files <outputs>`: ``<prefix>.contigs.fasta``,
    ``<prefix>.unitigs.fasta``, and ``<prefix>.unassembled.fasta``, where contigs are the primary
    output, unitigs are the primary output split at alternate paths,
    and unassembled are the leftover pieces.

    The :ref:`contigFilter <contigFilter>` parameter sets several parameters that control how small
    or low coverage initial contigs are handled.  By default, initial contigs with more than 50% of
    the length at less than 3X coverage will be classified as 'unassembled' and removed from the
    assembly, that is, ``contigFilter="2 0 1.0 0.5 3"``.  The filtering can be disabled by changing
    the last number from '3' to '0' (meaning, filter if 50% of the contig is less than 0X coverage).


Why is my assembly is missing my favorite short plasmid?
-------------------------------------
    In Canu v1.6 and earlier only the longest 40X of data (based on the specified genome size) is
    used for correction.  Datasets with uneven coverage or small plasmids can fail to generate
    enough corrected reads to give enough coverage for assembly, resulting in gaps in the genome or
    even no reads for small plasmids.  Set ``corOutCoverage=1000`` (or any value greater than your
    total input coverage) to correct all input data.

    An alternate approach is to correct all reads (``-correct corOutCoverage=1000``) then assemble
    40X of reads picked at random from the ``<prefix>.correctedReads.fasta.gz`` output.

    More recent Canu versions dynamically select poorly represented sequences to avoid missing short
    plasmids so this should no longer happen.

Why do I get less corrected read data than I asked for?
-------------------------------------
    Some reads are trimmed during correction due to being chimeric or because there wasn't enough
    evidence to generate a quality corrected sequence.  Typically, this results in a 25% loss.
    Setting ``corMinCoverage=0`` will report all bases, even low those of low quality.  Canu will
    trim these in its 'trimming' phase before assembly.


What is the minimum coverage required to run Canu?
-------------------------------------
    For HiFi data, we recommend 20-25x or higher. 
    For raw data, coverage more than 20X is typically enough to outperform current hybrid
    methods.  Below that, you will likely not assemble the full genome.  The following
    two papers have several examples.
     * `Koren et al. (2013) Reducing assembly complexity of microbial genomes with single-molecule sequencing <https://www.ncbi.nlm.nih.gov/pubmed/24034426>`_
     * `Koren and Walenz et al. (2017) Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation <https://www.ncbi.nlm.nih.gov/pubmed/28298431>`_

Can I use Illumina data too?
-------------------------------------
    No.  We've seen that using short reads for correction will homogenize repeats and
    mix up haplotypes.  Even though the short reads are very high quality, their length
    isn't sufficient for the true alignment to be identified, and so reads from other repeat
    instances are used for correction, resulting in incorrect corrections.

My circular element is duplicated/has overlap?
-------------------------------------
    This is expected for any circular elements. They can overlap by up to a read length due to how
    Canu constructs contigs. Canu provides an alignment string in the GFA output which can be
    converted to an alignment to identify the trimming points.

    An alternative is to run MUMmer to get self-alignments on the contig and use those trim
    points. For example, assuming the circular element is in ``tig00000099.fa``. Run::

      nucmer -maxmatch -nosimplify tig00000099.fa tig00000099.fa
      show-coords -lrcTH out.delta

    to find the end overlaps in the tig. The output would be something like::

      1	1895	48502	50400	1895	1899	99.37	50400	50400	3.76	3.77	tig00000001	tig00000001
      48502	50400	1	1895	1899	1895	99.37	50400	50400	3.77	3.76	tig00000001	tig00000001

    means trim to 1 to 48502. There is also an alternate `writeup
    <https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Circularizing-and-trimming>`_.

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

   Here is a quick walk-through using a command-line ftp client (should be available on most Linux
   and OSX installations). Say we want to transfer a file named ``reads.fastq``. First, run ``ftp
   ftp.cbcb.umd.edu``, specify ``anonymous`` as the user name and hit return for password
   (blank). Then ``cd incoming/sergek``, ``put reads.fastq``, and ``quit``.

   That's it, you won't be able to see the file but we can download it.
