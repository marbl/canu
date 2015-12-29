
.. _quickstart:

Canu Quick Start
================

Canu specializes in assembling PacBio or Oxford Nanopre sequences.  Canu will correct the reads, then trim suspicious regions (such as remaining SMRTbell adapter), then
assemble the corrected and cleaned reads into unitigs.

Brief Introduction
-------------------
Canu has been designed to auto-detect your resources and scale itself to fit. Two parameters let you restrict the resources used.

::

 maxMemory=XX
 maxThreads=XX

Memory is specified in gigabytes. On a single machine, it will restrict Canu to at most this limit, on the grid, no single job will try to use more than the specified resources.

The input sequences can be FASTA or FASTQ format, uncompressed, or compressed with gz, bz2 or xz.

Running on the grid
~~~~~~~~~~~~~~~~~~~~~~
Canu is designed to run on grid environments (LSF/PBS/Torque/Slrum/SGE are supported). Currently, Canu will submit itself to the default queue with default time options. You can overwrite this behavior by providing any specific parameters you want to be used for submission as an option. Users should also specify a job name to use on the grid:

::

 gridOptionsJobName=myassembly
 "gridOptions=--partition quick --time 2:00"

Assembling PacBio data
----------------------

Pacific Biosciences released P6-C4 chemistry reads.  You can download them
`directly <https://s3.amazonaws.com/files.pacb.com/datasets/secondary-analysis/e-coli-k12-P6C4/p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz>`_
(7 GB) or from the
`original page <https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly>`_.
You must have the Pac Bio SMRTpipe software installed to extract the reads as FASTQ.

We made a 25X subset FASTQ available
`here <http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq.gz>`_

or use the following curl command:

::

 curl -L -o p6.25x.fastq http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq
 
Correct, Trim and Assemble
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, canu will correct the reads, then trim the reads, then assemble the reads to unitigs.  

::

 canu \
  -p ecoli -d ecoli-auto \
  genomeSize=4.8m \
  -pacbio-raw p6.25x.fastq

This will use the prefix 'ecoli' to name files, compute the correction task in directory 'ecoli-auto/correction', the trimming task in directory 'ecoli-auto/trimming', and the unitig construction stage in 'ecoli-auto' itself.
Output files are described in the next section.

Find the Output
~~~~~~~~~~~~~~~~~~~~~~

Outputs from the assembly tasks are in:

- ecoli*/ecoli.layout
- ecoli*/ecoli.consensus.fasta

The canu progress chatter records statistics such as an input read histogram, corrected read histogram, and overlap types.


Correct, Trim and Assemble, Manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, however, it makes sense to do the three top-level tasks by hand.  This would allow trying
multiple unitig construction parameters on the same set of corrected and trimmed reads.

First, correct the raw reads::

 canu -correct \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   -pacbio-raw  p6.25x.fastq

Then, trim the output of the correction::

 canu -trim \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   -pacbio-corrected ecoli/correction/ecoli.correctedReads.fastq

And finally, assemble the output of trimming, twice::

 canu -assemble \
   -p ecoli -d ecoli-erate-0.025 \
   genomeSize=4.8m \
   errorRate=0.025 \
   -pacbio-corrected ecoli/trimming/ecoli.trimmedReads.fastq

 canu -assemble \
   -p ecoli -d ecoli-erate-0.035 \
   genomeSize=4.8m \
   errorRate=0.035 \
   -pacbio-corrected ecoli/trimming/ecoli.trimmedReads.fastq

The directory layout for correction and trimming is exactly the same as when we ran all tasks in the same command.
Each unitig construction task needs its own private work space, and in there the 'correction' and 'trimming' directories are empty.

Assembling Oxford Nanopore data
--------------------------------
A set of E. coli runs were released by the Loman lab.  You can download one
`directly <http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta>`_
or any of them from the
`original page <http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/>`_.

or use the following curl command:

::

 curl -L -o oxford.fasta http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta

Canu assembles any of the four available datasets into a single contig but we picked one dataset to use in this tutorial. Then, assemble the data as before::

 canu \
  -p ecoli -d ecoli-oxford \
  genomeSize=4.8m \
  -nanopore-raw oxford.fasta

The assembled identity is >99% before polishing.

Assembling With Multiple Technologies/Files 
-------------------------------------------

Canu takes an arbitrary number of input files/formats. We made a mixed dataset of about 10X of a PacBio P6 and 10X of an Oxford Nanopore run available `here <http://gembox.cbcb.umd.edu/mhap/raw/ecoliP6Oxford.tar.gz>`_

or use the following curl command:

::

 curl -L -o mix.tar.gz http://gembox.cbcb.umd.edu/mhap/raw/ecoliP6Oxford.tar.gz
 tar xvzf mix.tar.gz
 
Now you can assemble all the data::

 canu \
  -p ecoli -d ecoli-mix \
  genomeSize=4.8m \
  -pacbio-raw pacbio*fastq.gz \
  -nanopore-raw oxford.fasta.gz


Assembling Low Coverage Datasets
----------------------------------
When you have 30X or less coverage, it helps to adjust the Canu assembly parameters. You can download a 20X subset of `S. cerevisae <http://gembox.cbcb.umd.edu/mhap/raw/yeast_filtered.20x.fastq.gz>`_
 
or use the following curl command:

::

 curl -L -o yeast.20x.fastq.gz http://gembox.cbcb.umd.edu/mhap/raw/yeast_filtered.20x.fastq.gz

and run the assembler with modified parameters::

 canu \
  -p asm -d yeast \
  genomeSize=12.1m \
  corMhapSensitivity=high corMinCoverage=2 errorRate=0.035 \
  -pacbio-raw yeast.20x.fastq.gz
  

After the run completes, we can check the assembly statistics::

 tgStoreDump -sizes -s 12100000 -T yeast/unitigging/asm.tigStore 2 -G yeast/unitigging/asm.gkpStore

::

 lenSingleton n10 siz       7013 sum    1210884 idx        116
 lenSingleton sum    2338725 (genomeSize 12100000)
 lenSingleton num        416
 lenSingleton ave       5621
 lenAssembled n10 siz     696203 sum    1453015 idx          1
 lenAssembled n20 siz     575091 sum    2646269 idx          3
 lenAssembled n30 siz     550579 sum    3755422 idx          5
 lenAssembled n40 siz     455083 sum    5250476 idx          8
 lenAssembled n50 siz     392191 sum    6088423 idx         10
 lenAssembled n60 siz     205069 sum    7342769 idx         15
 lenAssembled n70 siz     140204 sum    8504891 idx         22
 lenAssembled n80 siz      99777 sum    9693133 idx         32
 lenAssembled n90 siz      64744 sum   10949303 idx         48
 lenAssembled n100 siz      15639 sum   12100894 idx         89
 lenAssembled sum   12607682 (genomeSize 12100000)
 lenAssembled num        150
 lenAssembled ave      84051

Consensus Accuracy
-------------------
While Canu corrects sequences and has 99% identity or greater with PacBio or Nanopore sequences, for the best accuracy we recommend polishing with a sequence-specific tool. We recommend `Quiver <http://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/HowToQuiver.rst>`_ for PacBio and `Nanopolish <http://github.com/jts/nanopolish>`_ for Oxford Nanpore data.

If you have Illumina sequences available, `Pilon <http://www.broadinstitute.org/software/pilon/>`_ can also be used to polish either PacBio or Oxford Nanopore assemblies.

Known Issues
-------------------

- LSF support has limited testing
- Large memory usage while unitig consensus calling on unitigs over 100MB in size (140Mb contig uses approximate 75GB).
- Distributed file systems (such as GPFS) causes issues with memory mapped files, slowing down parts of Canu, including meryl (0-mercounts) and falcon-sense (2-correction).
