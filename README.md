# Canu

Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page), designed for high-noise single-molecule sequencing (such as the [PacBio](http://www.pacb.com) [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/) or [Oxford Nanopore](https://www.nanoporetech.com/) [MinION](https://www.nanoporetech.com/products-services/minion-mki)).

Canu is a hierarchical assembly pipeline which runs in four steps:

* Detect overlaps in high-noise sequences using [MHAP](https://github.com/marbl/MHAP)
* Generate corrected sequence consensus
* Trim corrected sequences
* Assemble trimmed corrected sequences

## Install:

The easiest way to get started is to download a [release](http://github.com/marbl/canu/releases). 

Alternatively, you can also build the latest unreleased from github:

    git clone https://github.com/marbl/canu.git
    cd canu/src
    make -j <number of threads>

## Learn:

The [quick start](http://canu.readthedocs.io/en/stable/quick-start.html) will get you assembling quickly, while the [tutorial](http://canu.readthedocs.io/en/stable/tutorial.html) explains things in more detail.

## Run:

Brief command line help:

    ../<achitechture>/bin/canu

Full list of parameters:

    ../<architecture>/bin/canu -options

## Citation:
 - Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM. [Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation](http://dx.doi.org/10.1101/071282). bioRxiv. (2016).
