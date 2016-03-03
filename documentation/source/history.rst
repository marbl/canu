.. _history:

Software Background
====================
Canu is derived from the Celera Assembler. The Celera assembler [Myers 2000] 
was designed to reconstruct mammalian 
chromosomal DNA sequences from the short fragments of a whole genome shotgun 
sequencing project. The Celera Assembler was used to produce reconstructions of 
several large genomes, namely those of Homo sapiens [Venter 2001], Mus musculus 
[Mural 2002], Rattus norvegicus [unpublished data], Canis familiaris [Kirkness 
2003], Drosophila melanogaster [Adams 2000], and Anopheles gambiae [Holt 2001]. 
The Celera Assembler was shown to be very accurate when its reconstruction of 
the human genome was compared to independent reconstructions completed later 
[Istrail 2004]. It was used to reconstructing one of the first large-scale
metagenomic projects [Venter 2004, Rusch 2007] and a diploid human reference [Levy 2007, Denisov 2008]. 
It was adapted to 454 Pyrosequencing [Miller 2008] and PacBio sequencing [Koren 2011], demonstrating 
finished bacterial genomes [Koren 2013] and efficient algorithms for eukaryotic assembly [Berlin 2015].

In 2015 Canu was forked from Celera Assembler and specialized for single-molecule
high-noise sequences. The Celera Assembler codebase is no longer maintained.

Canu is a pipeline consisting of several executable programs and perl driver scripts.
The source code includes programs in the C++ language with Unix make scripts. The original
Celera Assembler was designed to run under Compaq(R) Tru64(R) Unix with access to 32GB 
RAM. It has also run under IBM(R) AIX(R) and Red Hat Linux. 

The Celera Assembler was released under the GNU General Public License, version 2 as as supplement
to the publication [Istrail 2004]. For the most recent license information please see
README.licences

References
--------------------
- Adams et al. (2000) The Genome Sequence of Drosophila melanogaster. Science 287 2185-2195.
- Myers et al. (2000) A Whole-Genome Assembly of Drosophila. Science 287 2196-2204.
- Venter et al. (2001) The Sequence of the Human Genome. Science 291 1304-1351.
- Mural et al. (2002) A Comparison of Whole-Genome Shotgun-Derived Mouse Chromosome 16 and the Human Genome. Science 296 1661-1671.
- Holt et al. (2002) The Genome Sequence of the Malaria Mosquito Anophelesd gambiae. Science 298 129-149.
- Istrail et al. (2004) Whole Genome Shotgun Assembly and Comparison of Human Genome Assemblies. PNAS 101 1916-1921.
- Kirkness et al. (2003) The Dog Genome: Survey Sequencing and Comparative Analysis. Science 301 1898-1903.
- Venter et al. (2004) Environmental genome shotgun sequencing of the Sargasso Sea. Science 304 66-74.
- Levy et al. (2007) The Diploid Genome Sequence of an Individual Human. PLoS Biology 0050254
- Rusch et al. (2007) The Sorcerer II Global Ocean Sampling Expedition: Northwest Atlantic through Eastern Tropical Pacific. PLoS Biology 1821060.
- Denisov et al. (2008) Consensus Generation and Variant Detection by Celera Assembler. Bioinformatics 24(8):1035-40
- Miller et al. (2008) Aggressive Assembly of Pyrosequencing Reads with Mates. Bioinformatics 24(24):2818-2824
- Koren et al. (2012) Hybrid error correction and de novo assembly of single-molecule sequencing reads, Nature Biotechnology, July 2012.
- Koren et al. (2013) Reducing assembly complexity of microbial genomes with single-molecule sequencing, Genome Biology 14:R101.
- Berlin et. al. (2015) Assembling Large Genomes with Single-Molecule Sequencing and Locality Sensitive Hashing. Nature Biotechnology. (2015).
