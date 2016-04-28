# Canu

Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page), designed for high-noise single-molecule sequencing (such as the [PacBio](http://www.pacb.com) [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/) or [Oxford Nanopore](https://www.nanoporetech.com/) [MinION](https://www.nanoporetech.com/products-services/minion-mki)).

Canu is a hierachical assembly pipeline which runs in four steps:

* Detect overlaps in high-noise sequences using [MHAP](https://github.com/marbl/MHAP)
* Generate corrected sequence consensus
* Trim corrected sequences
* Assemble trimmed corrected sequences

## Build:

    git clone https://github.com/marbl/canu.git
    cd canu/src
    make -j <number of threads>

## Run:

Brief command line help:

    ../<achitechture>/bin/canu
    

Full list of parameters:

    ../<architecture>/bin/canu -options
    
## Learn:

The [quick start](http://canu.readthedocs.io/en/stable/quick-start.html) will get you assembling quickly, while the [tutorial](http://canu.readthedocs.io/en/stable/tutorial.html) explains things in more detail.

## Citation:

 - Berlin K, Koren S, Chin CS, Drake PJ, Landolin JM, Phillippy AM [Assembling Large Genomes with Single-Molecule Sequencing and Locality Sensitive Hashing](http://www.nature.com/nbt/journal/v33/n6/abs/nbt.3238.html). Nature Biotechnology. (2015).
 - Stay tuned for a Canu-specific citation
