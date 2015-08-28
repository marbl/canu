
.. _quickstart:

Canu Quick Start
================

Canu specializes in assembling PacBio filtered_subreads directly, without correction.  To be fair,
canu will correct the reads, then trim suspicious regions (such as remaining SMRTbell adapter), then
assemble the corrected and cleaned reads into unitigs.

Download the Reads
-------------------

There are two datasets for Escherichia coli K12 MG1655.

Our friends at NBACC released the e.coli reads used in
`Koren et al. <http://genomebiology.com/2013/14/9/R101>`_ as NCBI SRA project
`SRP020003 <[http://www.ncbi.nlm.nih.gov/sra/?term=SRP020003>`_.  The SRA contains reads generated using
`454 FLX Titanium <http://www.ncbi.nlm.nih.gov/sra/SRX255226>`_,
`PacBio RS C2 CCS <http://www.ncbi.nlm.nih.gov/sra/SRX255779>`_ and
`PacBio RS C2 <http://www.ncbi.nlm.nih.gov/sra/SRX255228>`_.
The full data set is 600x coverage (11.4 GB), and the reads are (slightly) older technology, so
we'll not use them in this example. If we were to use them, we could filter out all the short reads,
or sample 100x randomly from all files, or just use one of the libraries.

Pacific Biosciences released 100x of P4-C2 chemistry reads.  You can download them
`directly <http://files.pacb.com/datasets/secondary-analysis/ecoli-k12-P4C2-20KSS/ecoliK12.tar.gz>`_
(6.35 GB) or from the
`original page <https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-20kb-Size-Selected-Library-with-P4-C2>`_.
You must have the Pac Bio SMRTpipe software installed to extract the reads as FASTQ.

We made FASTQ files for
`Celera Assembler <http://wgs-assembler.sourceforge.net/wiki/index.php/Escherichia_coli_K12_MG1655,_using_uncorrected_PacBio_reads,_with_CA8.2>`_
that are here:

- `escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz <http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz/download>`_ (116 MB)
- `escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.2.fastq.xz <http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.2.fastq.xz/download>`_ (110 MB)
- `escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.3.fastq.xz <http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.3.fastq.xz/download>`_ (136 MB))

or use the following impressively long curl commands:

::

 curl -L -o escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz/download
 curl -L -o escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.2.fastq.xz http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.2.fastq.xz/download
 curl -L -o escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.3.fastq.xz http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.0/datasets/escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.3.fastq.xz/download

Correct, Trim and Assemble
--------------------------

By default, canu will correct the reads, then trim the reads, then assemble the reads to unitigs.  Since this is a quick start, we'll use only one of the three files.

::

 canu \
  -p ecoli -d ecoli-auto \
  genomeSize=4.8m \
  errorRate=0.02 \
  -pacbio-raw escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz

This will use the prefix 'ecoli' to name files, compute the correction task in directory 'ecoli-auto/correction', the trimming task in directory 'ecoli-auto/trimming', and the unitig construction stage in 'ecoli-auto' itself.
Output files are described in the next section.

Correct, Trim and Assemble, Manually
------------------------------------

Sometimes, however, it makes sense to do the three top-level tasks by hand.  This would allow trying
multiple unitig construction parameters on the same set of corrected and trimmed reads.

First, correct the raw reads::

 canu -correct \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   errorRate=0.02 \
   -pacbio-raw  escherichia_coli_k12_mg1655.m130404_014004_sidney_c100506902550000001823076808221337_s1_p0.1.fastq.xz

Then, trim the output of the correction::

 canu -trim \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   errorRate=0.02 \
   -pacbio-corrected ecoli/correction/ecoli.correctedReads.fastq

And finally, assemble the output of trimming, twice::

 canu -assemble \
   -p ecoli -d ecoli-erate-0.02 \
   genomeSize=4.8m \
   errorRate=0.02 \
   -pacbio-corrected ecoli/trimming/ecoli.trimmedReads.fastq

 canu -assemble \
   -p ecoli -d ecoli-erate-0.03 \
   genomeSize=4.8m \
   errorRate=0.03 \
   -pacbio-corrected ecoli/trimming/ecoli.trimmedReads.fastq

The directory layout for correction and trimming is exactly the same as when we ran all tasks in the same command.
Each unitig construction task needs its own private work space, and in there the 'correction' and 'trimming' directories are empty.

Find the Output
---------------

Outputs from the assembly tasks are in:

- ecoli*/ecoli.layout
- ecoli*/ecoli.consensus.fastq

though canu doesn't report this in a nice, easy to find way, yet -- the output files are buried in the canu progress chatter.  Nor are statistics output.


