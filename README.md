# Canu

Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page), designed for high-noise single-molecule sequencing (such as the [PacBio](http://www.pacb.com) [RS II](http://www.pacb.com/products-and-services/pacbio-systems/rsii/)/[Sequel](http://www.pacb.com/products-and-services/pacbio-systems/sequel/) or [Oxford Nanopore](https://www.nanoporetech.com/) [MinION](https://nanoporetech.com/products)).

Canu is a hierarchical assembly pipeline which runs in four steps:

* Detect overlaps in high-noise sequences using [MHAP](https://github.com/marbl/MHAP)
* Generate corrected sequence consensus
* Trim corrected sequences
* Assemble trimmed corrected sequences

## Install:

* Do NOT download the .zip source code.  It is missing files and will not compile.  This is a [known flaw](https://github.com/dear-github/dear-github/issues/214) with git itself.

* The easiest way to get started is to download a binary [release](http://github.com/marbl/canu/releases).

* Installing with a 'package manager' is not encouraged, but if you have no other choice:
  * Conda: `conda install -c conda-forge -c bioconda -c defaults canu`
  * Homebrew: `brew install brewsci/bio/canu`

* Alternatively, you can use the latest unreleased version from the source
  code.  This version has not undergone the same testing as a release and so
  may have unknown bugs or issues generating sub-optimal assemblies. We
  recommend the release version for most users.

        git clone https://github.com/marbl/canu.git
        cd canu/src
        make -j <number of threads>

  * FreeBSD generally requires libboost be installed from packages/ports.  It
    will compile with either clang (>= 14) or gcc (>= 9).  It requires
    openjdk18.

        With clang, (default 14) needs libboost from ports.

          gmake

        With gcc (9+), can use the canu-supplied libboost or libboost from ports.

          gmake CC=gcc9 CXX=g++9 BOOST=libboost    #  Canu-supplied boost
          gmake CC=gcc9 CXX=g++9                   #  Ports/packages supplied boost

  * MacOS Apple Silicon requires libboost, and either openjdk or oracle-jdk
    to be installed from homebrew (preferred) or MacPorts.  It will compile
    with either clang (>=14) or gcc (>= 9) but WILL NOT compile with the
    standard Xcode compiler.

          make CC=gcc-11 CXX=g++-11 BOOST=libboost  #  Ports/packages supplied boost
          make CC=gcc-11 CXX=g++-11                 #  Ports/packages supplied boost

  * MacOS Intel is probably the same as Apple Silicon, but not tested.

  * Linux does not need a system installed libboost.

* An *unsupported* Docker image made by Frank Förster is at https://hub.docker.com/r/greatfireball/canu/.

## Learn:

The [quick start](http://canu.readthedocs.io/en/latest/quick-start.html) will get you assembling quickly, while the [tutorial](http://canu.readthedocs.io/en/latest/tutorial.html) explains things in more detail.

## Run:

Brief command line help:

    ../<architecture>/bin/canu

Full list of parameters:

    ../<architecture>/bin/canu -options

## Citation:
 - Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM. [Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation](https://doi.org/10.1101/gr.215087.116). Genome Research. (2017). `doi:10.1101/gr.215087.116`
 - Koren S, Rhie A, Walenz BP, Dilthey AT, Bickhart DM, Kingan SB, Hiendleder S, Williams JL, Smith TPL, Phillippy AM. [De novo assembly of haplotype-resolved genomes with trio binning](http://doi.org/10.1038/nbt.4277).  Nature Biotechnology.  (2018). (If you use trio-binning)
 - Nurk S, Walenz BP, Rhiea A, Vollger MR, Logsdon GA, Grothe R, Miga KH, Eichler EE, Phillippy AM, Koren S. [HiCanu: accurate assembly of segmental duplications, satellites, and allelic variants from high-fidelity long reads](https://doi.org/10.1101/2020.03.14.992248).  biorXiv.  (2020). (If you use -pacbio-hifi)
 

