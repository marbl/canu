# Canu

<p><img style="float: right; margin: 0px 15px 15px 0px;" src="https://raw.githubusercontent.com/marbl/canu/master/logo.jpg" width="100" /> Canu is a fork of the [Celera Assembler](http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page "Celera Assembler") designed for high-noise single-molecule sequencing (such as the PacBio RSII or Oxford Nanopore MinION). The software is currently alpha level, feel free to use and report issues encountered.
<br style="clear: both;" />
</p>

## Build

    git clone https://github.com/marbl/canu.git
    cd canu/src
    make
    
For a quick user-quide, run:

    ../<achitechture>/bin/canu
    

For full list of options, run:
    ../<architecture>/bin/canu -options
    
## Docs
Canu is a hierachical assembly pipeline which runs in three steps:

* Detect overlap in high-noise sequences using [MHAP](https://github.com/marbl/MHAP "MHAP")
* Generate corrected sequence consensus
* Assemble corrected sequences

Stay tuned, more coming soon.
