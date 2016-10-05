
.. _faq:

Canu FAQ
========================

**Q**:
    What resources does Canu require for a bacterial genome assembly? A mammalian assembly?

**A**:
    Canu is designed to scale resources to the system it runs on. It will report if the a system does not meet the minimum requirements for a given genome size.
    
    Typically, a bacterial genome can be assembled in 1-10 cpu hours, depending on coverage (~20 min on 16-cores) and 4GB of ram (8GB is recommended). A mammalian genome (such as human) can be assembled in 10-25K cpu hours, depending on coverage (a grid environment is recommended) and at least one machine with 64GB of ram (128GB is recommended).
    
**Q**:
    What parameters should I use for my genome? Sequencing type?
    
**A**:
    By default, Canu is designed to be universal on a large range of PacBio (C2-P6-C4) and Oxford Nanopore (R6-R9) data. You can adjust parameters to increase efficiency for your datatype:
    
    **Nanopore R7 1D**
       With R7 1D sequencing data, multiple rounds of error correction are helpful. Generally, this is true if you input data has lower identity than 80%. You can run just the correction from Canu with the options
    
    ::
    
        -correct corOutCoverage=500 corMinCoverage=0 corMhapSensitivity=high
    
    for 5-10 rounds, supplying the asm.correctedReads.fasta.gz output from round ``i-1`` to round ``i``. Assemble with
    
    ::
    
        -nanopore-corrected <your data> errorRate=0.1 utgGraphDeviation=50    
     **Nanopore R7 2D**
        ``errorRate=0.025``
    **Nanopore R9 1D**
        ``errorRate=0.025``
    **Nanopore R9 2D**
        ``errorRate=0.013``
    **PacBio P6**
        ``errorRate=0.013``
    **PacBio Sequel** (based on *A. thaliana* `publically released  dataset <http://www.pacb.com/blog/sequel-system-data-release-arabidopsis-dataset-genome-assembly/>`_)
        ``errorRate=0.013 corMhapSensitivity=normal``

**Q**:
    My assembly continuity is not good, how can I improve it?
    
**A**:
    The most important determinant for assembly quality is sequence length, followed by the repeat complexity/heterozygosity of your sample.  The first thing to check is the amount of corrected bases output by the correction step.  This is logged in the stdout of Canu or in canu-scripts/canu.*.out if you are running in a grid environment. For example on `a haploid H. sapiens <https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN02744161>`_ sample:
    
    ::
    
       -- BEGIN TRIMMING
       --
       ...
       -- In gatekeeper store 'chm1/trimming/asm.gkpStore':
       --   Found 5459105 reads.
       --   Found 91697412754 bases (29.57 times coverage).
       --
       ...

   Canu tries to correct the longest 40X of data. Some loss is normal but having output coverage below 20-25X is a sign that correction did not work well (assuming you have more input coverage than that). If that is the case, re-running with ``corMhapSensitivity=normal`` if you have >50X or ``corMhapSensitivity=high corMinCoverage=0`` otherwise can help. You can also increase the target coverage to correct ``corOutCoverage=100`` to get more correct sequences for assembly. If there are sufficient corrected reads, the poor assembly is likely due to either  repeats in the genome being greater than read lengths or a high heterozygosity in the sample. Stay tuned for mor information on tuning unitigging in those instances.

**Q**:
    How do I run Canu on my SLURM/SGE/PBS/LSF/Torque system?

**A**:
    Canu will auto-detect and configure itself to submit on most grids. If your grid requires special options (such as a partition on SLURM or an account code on SGE, specify it with ``gridOptions="<your options list>"`` which will passed to the sheduler by Canu. If you have a grid system but prefer to run locally, specify useGrid=false
    
**Q**:
    My asm.contigs.fasta is empty, why?

**A**:
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

    The default filtering is ``2 1000 0.75 0.75 2``. If you are assembling amplified data or viral data, it is possible your assembly will be flagged as unassembled. In those cases, you can turn off the filtering with the parameters

    ::

       contigFilter="2 1000 1.0 1.0 2"

**Q**:
    Why is my assembly is missing my favorite short plasmid X?

**A**:
    The first step in Canu is to find high-error overlaps and generate corrected sequences for subsequent assembly. This is currently the fastest step in Canu. By default, only the longest 40X of data (based on the specified genome size) is used for correction. If you have a dataset with uneven coverage or small plasmids, correcting the longest 40X may not give you sufficient coverage of your genome/plasmid. In these cases, you can set

    ::

        corOutCoverage=1000

    Or any large value greater than your total input coverage which will correct and assemble all input data, at the expense of runtime. This option is also recommended for metagenomic datasets where all data is useful for assembly.

**Q**:
    Why do I get only 30X of corrected data?

**A**:
    By default, only the longest 40X of data (based on the specified genome size) is used for correction. Typically, some reads are trimmed during correction due to being chimeric or having erroneous sequence, resulting in a loss of 20-25% (30X output). You can force correction to be non-lossy by setting

    ::  

       corMinCoverage=0

    In which case the corrected reads output will be the same length as the input data, keeping any high-error unsupported bases. Canu will trim these in downstream steps before assembly.

**Q**:
   What is the minimum coverage required to run Canu?

**A**:
    We have found that on eukaryotic genomes >=20X typically begins to outperform current hybrid methods. For low coverage datasets (<=30X) we recommend the following parameters

    ::

       corMinCoverage=0 errorRate=0.035

    For high-coverage datasets (typically >=60X) you can decrease the error rate since the higher number of reads should allow sufficient assembly from only the best subset

    ::

       errorRate=0.013

    However, the above is mainly an optimization for speed and will not affect your assembly continuity.


**Q**:
   My genome is AT/GC rich, do I need to adjust parameters?

**A**:
    On bacterial genomes, typically no. On repetitive genomes with AT<=25 or 75>=AT (or GC) the sequence biases the Jaccard estimate used by MHAP. In those cases setting

    ::

        corMaxEvidenceErate=0.15

    has been sufficient to correct for the bias in our testing. In general, with high coverage repetitive genomes (such as plants) it can be beneficial to set the above parameter as it will eliminate repetitive matches, speed up the assembly, and sometime improve unitigs.
