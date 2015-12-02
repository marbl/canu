
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

Running on the grid
~~~~~~~~~~~~~~~~~~~~~~
Canu is designed to run on grid environments (LSF/PBS/Torque/Slrum/SGE are supported). Currently, Canu will submit itself to the default queue with default time options. You can overwrite this behavior by providing any specific parameters you want to be used for submission as an option. Users should also specify a job name to use on the grid:

::

 gridOptionsJobName=myassembly
 "gridOptions=--partition quick --time 2:00"

Assembling PacBio data
-------------------

Pacific Biosciences released P6-C4 chemistry reads.  You can download them
`directly <https://s3.amazonaws.com/files.pacb.com/datasets/secondary-analysis/e-coli-k12-P6C4/p6c4_ecoli_RSII_DDR2_with_15kb_cut_E01_1.tar.gz>`_
(7 GB) or from the
`original page <https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly>`_.
You must have the Pac Bio SMRTpipe software installed to extract the reads as FASTQ.

We made a 25X subset FASTQ available
`here <http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq.gz>`_

or use the following curl command:

::

 curl -L -o p6.25x.fastq.gz http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq.gz
 
Correct, Trim and Assemble
~~~~~~~~~~~~~~~~~~~~~~

By default, canu will correct the reads, then trim the reads, then assemble the reads to unitigs.  

::

 canu \
  -p ecoli -d ecoli-auto \
  genomeSize=4.8m \
  -pacbio-raw p6.25x.fastq.gz

This will use the prefix 'ecoli' to name files, compute the correction task in directory 'ecoli-auto/correction', the trimming task in directory 'ecoli-auto/trimming', and the unitig construction stage in 'ecoli-auto' itself.
Output files are described in the next section.

Correct, Trim and Assemble, Manually
~~~~~~~~~~~~~~~~~~~~~~

Sometimes, however, it makes sense to do the three top-level tasks by hand.  This would allow trying
multiple unitig construction parameters on the same set of corrected and trimmed reads.

First, correct the raw reads::

 canu -correct \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   -pacbio-raw  p6.25x.fastq.gz

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

Find the Output
~~~~~~~~~~~~~~~~~~~~~~

Outputs from the assembly tasks are in:

- ecoli*/ecoli.layout
- ecoli*/ecoli.consensus.fasta

The canu progress chatter records statistics such as an input read histogram, corrected read histogram, and overlap types.

Assembling Oxford Nanopore data
-------------------
A set of E. coli runs were released by the Loman lab.  You can download them
`directly <http://nanopore.s3.climb.ac.uk/MAP006-2_2D_pass.fasta>`_
(7 GB) or from the
`original page <http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/>`_.

or use the following curl command:

::

 curl -L -o oxford.fasta http://nanopore.s3.climb.ac.uk/MAP006-2_2D_pass.fasta

Canu assembles any of the four available datasets into a single contig but we picked one dataset to use in this tutorial. Then, assemble the data as before::

 canu \
  -p ecoli -d ecoli-oxford \
  genomeSize=4.8m \
  -nanopore-raw oxford.fasta

Assembling Low Coverage
-------------------
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

 tgStoreDump -stats -s 12100000 -T yeast/unitigging/asm.tigStore 2 -G yeast/unitigging/asm.gkpStore

Known Issues
-------------------

- LSF support has limited testing
- Large virtual memory usage while generating corrected sequences
- Large memory usage while unitig consensus calling on contigs over 50MB in size
