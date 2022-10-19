
.. _quickstart:

Canu Quick Start
================

Canu specializes in assembling PacBio or Oxford Nanopore sequences.  Canu operates in three phases:
correction, trimming and assembly.  The correction phase will improve the accuracy of bases in
reads.  The trimming phase will trim reads to the portion that appears to be high-quality sequence,
removing suspicious regions such as remaining SMRTbell adapter.  The assembly phase will order the
reads into contigs, generate consensus sequences and create graphs of alternate paths.

For eukaryotic genomes, coverage more than 20x is enough to outperform current hybrid methods,
however, between 30x and 60x coverage is the recommended minimum.  More coverage will let Canu use
longer reads for assembly, which will result in better assemblies.

Input sequences can be FASTA or FASTQ format, uncompressed or compressed with gzip (.gz), bzip2
(.bz2) or xz (.xz).  Note that zip files (.zip) are not supported.

Canu can resume incomplete assemblies, allowing for recovery from system outages or other abnormal
terminations.  On each restart of Canu, it will examine the files in the assembly directory to
decide what to do next.  For example, if all but two overlap tasks have finished, only the two that
are missing will be computed.  For best results, do not change Canu parameters between restarts.

Canu will auto-detect computational resources and scale itself to fit, using all of the resources
available and are reasonable for the size of your assembly.  Memory and processors can be explicitly
limited with with parameters :ref:`maxMemory <maxMemory>` and :ref:`maxThreads <maxThreads>`.  See section :ref:`execution`
for more details.

Canu will automatically take full advantage of any LSF/PBS/PBSPro/Torque/Slrum/SGE grid available,
even submitting itself for execution.  Canu makes heavy use of array jobs and requires job
submission from compute nodes, which are sometimes not available or allowed.  Canu option
``useGrid=false`` will restrict Canu to using only the current machine, while option
``useGrid=remote`` will configure Canu for grid execution but not submit jobs to the grid.
See section :ref:`execution` for more details.

The :ref:`tutorial` has more background, and the :ref:`faq` has a wealth of practical advice.


Assembling PacBio CLR or Nanopore data
----------------------

Pacific Biosciences released P6-C4 chemistry reads for Escherichia coli K12.  You can `download them
from their original release
<https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly>`_, but note that you
must have the `SMRTpipe software <http://www.pacb.com/support/software-downloads/>`_ installed to
extract the reads as FASTQ.  Instead, use a `FASTQ format 25X subset
<http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq>`_ (223MB).  Download from the command line
with::

 curl -L -o pacbio.fastq http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq

Confirm the MD5SUM matches::

 9bb4c10c41c5442d630af8b504042334  pacbio.fastq

There doesn't appear to be any "official" Oxford Nanopore sample data, but the `Albertsen Lab <https://albertsenlab.org/>`_ released a `run
<https://albertsenlab.org/we-ar10-3-pretty-close-now/>`_, also for Escherichia coli K12.  Download the R10 data from `FigShare <https://figshare.com/articles/dataset/Ecoli_K12_MG1655_R10_3_HAC/11823087>`_

Confirm the MD5SUM matches::

 14ae7c31805ab048e0b413261857d82f ecolk12mg1655_R10_3_guppy_345_HAC.fastq.gz

By default, Canu will correct the reads, then trim the reads, then assemble the reads to unitigs.
Canu needs to know the approximate genome size (so it can determine coverage in the input reads)
and the technology used to generate the reads.

For PacBio::

 canu \
  -p ecoli -d ecoli-pacbio \
  genomeSize=4.8m \
  -pacbio pacbio.fastq

For Nanopore::

 canu \
  -p ecoli -d ecoli-oxford \
  genomeSize=4.8m maxInputCoverage=100 \
  -nanopore ecolk12mg1655_R10_3_guppy_345_HAC.fastq

Output and intermediate files will be in directories 'ecoli-pacbio' and 'ecoli-nanopore',
respectively.  Intermediate files are written in directories 'correction', 'trimming' and
'unitigging' for the respective stages.  Output files are named using the '-p' prefix, such as
'ecoli.contigs.fasta', 'ecoli.unitigs.gfa', etc.  See section :ref:`outputs` for more details on
outputs (intermediate files aren't documented).

Assembling PacBio HiFi with HiCanu
----------------------

HiCanu has support for PacBio HiFi data by compressing homopolymers, correcting isolated errors, and masking systematic errors. We will now assemble and `E. coli K12
<https://sra-pub-src-1.s3.amazonaws.com/SRR10971019/m54316_180808_005743.fastq.1>`_ HiFi dataset sequenced by PacBio available at `NCBI SRA <https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10971019>`_ (3 GB).  When assembling, we
use `-pacbio-hifi` to specify the input reads::

 curl -L -o ecoli.fastq https://sra-pub-src-1.s3.amazonaws.com/SRR10971019/m54316_180808_005743.fastq.1

 canu \
  -p asm -d ecoli_hifi \
  genomeSize=4.8m \
  -pacbio-hifi ecoli.fastq
  
Trio Binning Assembly
----------------------------------

Canu has support for using parental short-read sequencing to classify and bin the F1 reads (see `Trio Binning manuscript
<https://www.biorxiv.org/content/early/2018/02/26/271486>`_ for details). This example demonstrates the functionality using a synthetic mix of two Escherichia coli datasets.  First download the data ::

 curl -L -o K12.parental.fasta https://gembox.cbcb.umd.edu/triobinning/example/k12.12.fasta
 curl -L -o O157.parental.fasta https://gembox.cbcb.umd.edu/triobinning/example/o157.12.fasta
 curl -L -o F1.fasta https://gembox.cbcb.umd.edu/triobinning/example/pacbio.fasta

Confirm the MD5SUM matches::

 69920456a2ef25fc3e89cdcb604861ed  K12.parental.fasta
 792d0af43740b3534516e8f73ead8a35  O157.parental.fasta
 64c8befea83d043344bdff7c43b04a71  F1.fasta

and run Canu::

 canu \
  -p asm -d ecoliTrio \
  genomeSize=5m \
  -haplotypeK12 K12.parental.fasta \
  -haplotypeO157 O157.parental.fasta \
  -pacbio F1.fasta

The run will first bin the reads into the haplotypes (``ecoliTrio/haplotype/haplotype-*.fasta.gz``) and provide a summary of the classification in ``ecoliTrio/haplotype/haplotype.log``::

  -- Processing reads in batches of 100 reads each.
  --
  --   119848 reads    378658103 bases written to haplotype file ./haplotype-K12.fasta.gz.
  --   308353 reads   1042955878 bases written to haplotype file ./haplotype-O157.fasta.gz.
  --     4114 reads      6520294 bases written to haplotype file ./haplotype-unknown.fasta.gz.


Next, the haplotypes are assembled in ``ecoliTrio/asm-haplotypeK12/asm-haplotypeK12.contigs.fasta`` and ``ecoliTrio/asm-haplotypeO157/asm-haplotypeO157.contigs.fasta``. By default, if the unassigned bases are > 5% of the total, they are included in both haplotypes. This can be controlled with the :ref:`hapUnknownFraction <hapUnknownFraction>` option. 

As comparison, you can try co-assembling the datasets instead::

 canu \
  -p asm -d ecoliHap \
  genomeSize=5m \
  corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \
 -pacbio F1.fasta

and compare the continuity/accuracy. 

Please note, trio binning is designed to work with raw sequences prior to correction. Do not correct the reads together and then run trio-binning, this will not work and Canu will give an error.

Trio binning does not yet support inputting PacBio HiFi reads for binning as they get flagged as "corrected" and the same error as above is given. As a workaround, run ``canu -haplotype`` specifying the HiFi reads as -pacbio-raw. This will bin the data and create shell scripts to start the assembly. Edit the shell scripts to replace -pacbio-raw with -pacbio-corrected or -pacbio-hifi and run the assemblies manually.

Assembling With Multiple Technologies and Multiple Files
-------------------------------------------

Canu can use reads from any number of input files, which can be a mix of formats and technologies. Note that current combining PacBio HiFi data with other datatypes it not supported. We'll assemble a mix of 10X PacBio CLR reads in two FASTQ files and 10X of Nanopore reads in one FASTA
file::

 curl -L -o mix.tar.gz http://gembox.cbcb.umd.edu/mhap/raw/ecoliP6Oxford.tar.gz
 tar xvzf mix.tar.gz

 canu \
  -p ecoli -d ecoli-mix \
  genomeSize=4.8m \
  -pacbio pacbio.part?.fastq.gz \
  -nanopore oxford.fasta.gz


Correct, Trim and Assemble, Manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, however, it makes sense to do the three top-level tasks by hand.  This would allow trying
multiple unitig construction parameters on the same set of corrected and trimmed reads, or skipping
trimming and assembly if you only want corrected reads.

We'll use the PacBio reads from above.  First, correct the raw reads::

 canu -correct \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   -pacbio  pacbio.fastq

Then, trim the output of the correction::

 canu -trim \
   -p ecoli -d ecoli \
   genomeSize=4.8m \
   -corrected -pacbio ecoli/ecoli.correctedReads.fasta.gz

And finally, assemble the output of trimming, twice, with different stringency on which overlaps to
use (see :ref:`correctedErrorRate <correctedErrorRate>`)::

 canu \
   -p ecoli -d ecoli-erate-0.039 \
   genomeSize=4.8m \
   correctedErrorRate=0.039 \
   -trimmed -corrected -pacbio ecoli/ecoli.trimmedReads.fasta.gz

 canu \
   -p ecoli -d ecoli-erate-0.075 \
   genomeSize=4.8m \
   correctedErrorRate=0.075 \
   -trimmed -corrected -pacbio ecoli/ecoli.trimmedReads.fasta.gz

Note that the assembly stages use different '-d' directories.  It is not possible to run multiple
copies of canu with the same work directory.

You can also try uncorrected ONT assembly which works for higher quality (95% accuracy) data, though this mode should be considered experimental::

 canu \
  -p ecoli -d ecoli-oxford-uncorrected \
  genomeSize=4.8m \
  -untrimmed correctedErrorRate=0.12 maxInputCoverage=100 'batOptions=-eg 0.10 -sb 0.01 -dg 2 -db 1 -dr 3' -pacbio-hifi ecolk12mg1655_R10_3_guppy_345_HAC.fastq 


Assembling Low Coverage Datasets
----------------------------------

We claimed Canu works down to 20X coverage, and we will now assemble `a 20X subset of S. cerevisae
<http://gembox.cbcb.umd.edu/mhap/raw/yeast_filtered.20x.fastq.gz>`_ (215 MB).  When assembling, we
adjust :ref:`correctedErrorRate <correctedErrorRate>` to accommodate the slightly lower
quality corrected reads::

 curl -L -o yeast.20x.fastq.gz http://gembox.cbcb.umd.edu/mhap/raw/yeast_filtered.20x.fastq.gz

 canu \
  -p asm -d yeast \
  genomeSize=12.1m \
  correctedErrorRate=0.105 \
  -pacbio yeast.20x.fastq.gz

Consensus Accuracy
-------------------

HiCanu consensus sequences using PacBio HiFi data are typically well above 99.99% We discourage any post-processing/polishing of these assemblies as mis-mapping within repeats can introduce errors.

Canu consensus sequences are typically well above 99% identity for PacBio datasets.  Nanopore accuracy varies depending on pore and basecaller version, but is typically above 99% for recent data. Accuracy can be improved by
polishing the contigs with tools developed specifically for that task.  We recommend `Arrow
<http://github.com/PacificBiosciences/GenomicConsensus>`_ for PacBio and `Nanopolish
<http://github.com/jts/nanopolish>`_ or `Medaka <https://github.com/nanoporetech/medaka>`_ for Oxford Nanpore data.
When Illumina reads are available, `FreeBayes <https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish>`_
can be used to polish either PacBio or Oxford Nanopore assemblies.
