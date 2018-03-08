.. _history:

Software Background
====================

Canu is derived from `Celera Assembler <http://wgs-assembler.sourceforge.net/>`_, which is no longer maintained.

Celera Assembler [`Myers 2000 <http://doi.org/10.1126/science.287.5461.2196>`_] was designed to reconstruct mammalian chromosomal DNA
sequences from the short fragments of a whole genome shotgun sequencing project.
Celera Assembler was used to produce reconstructions of several large genomes, namely
those of Homo sapiens [`Venter 2001 <http://doi.org/10.1126/science.1058040>`_], Mus musculus [`Mural 2002 <http://doi.org/10.1126/science.1069193>`_], Rattus norvegicus
[`unpublished data <https://www.genome.gov/10001855/rat-genome-sequencing/>`_], Canis familiaris [`Kirkness 2003 <http://doi.org/10.1126/science.1086432>`_], Drosophila melanogaster [`Adams
2000 <http://doi.org/10.1126/science.287.5461.2185>`_], and Anopheles gambiae [`Holt 2001 <http://doi.org/10.1126/science.1076181>`_].  Celera Assembler was shown to be very
accurate when its reconstruction of the human genome was compared to independent
reconstructions completed later [`Istrail 2004 <http://doi.org/10.1073/pnas.0307971100>`_]. It was used to reconstructing one of the
first large-scale metagenomic projects [`Venter 2004 <http://doi.org/10.1126/science.1093857>`_, `Rusch 2007 <http://doi.org/10.1371/journal.pbio.0050077>`_] and a diploid human
reference [`Levy 2007 <http://doi.org/10.1371/journal.pbio.0050254>`_, `Denisov 2008 <http://doi.org/10.1093/bioinformatics/btn074>`_].  It was adapted to 454 Pyrosequencing [`Miller 2008 <http://doi.org/10.1093/bioinformatics/btn548>`_]
and PacBio sequencing [`Koren 2012 <http://doi.org/10.1038/nbt.2280>`_], demonstrating finished bacterial genomes [`Koren 2013 <http://doi.org/10.1186/gb-2013-14-9-r101>`_]
and efficient algorithms for eukaryotic assembly [`Berlin 2015 <http://doi.org/10.1038/nbt.3238>`_].

Celera Assembler was released under the GNU General Public License, version 2 as a
supplement to [`Istrail 2004 <http://doi.org/10.1073/pnas.0307971100>`_].

Canu [`Koren and Walenz 2017 <http://doi.org/10.1101/gr.215087.116>`_] was branched from
Celera Assembler in 2015, and specialized for single-molecule high-noise sequences.
For the most recent license information on Canu,
see `README.licences <https://github.com/marbl/canu/blob/master/README.licenses>`_.

References
--------------------
- Adams et al.            (2000) `The Genome Sequence of Drosophila melanogaster                                                        <http://doi.org/10.1126/science.287.5461.2185>`_. Science 287 2185-2195.
- Myers et al.            (2000) `A Whole-Genome Assembly of Drosophila                                                                 <http://doi.org/10.1126/science.287.5461.2196>`_. Science 287 2196-2204.
- Venter et al.           (2001) `The Sequence of the Human Genome                                                                      <http://doi.org/10.1126/science.1058040>`_.       Science 291 1304-1351.
- Mural et al.            (2002) `A Comparison of Whole-Genome Shotgun-Derived Mouse Chromosome 16 and the Human Genome                 <http://doi.org/10.1126/science.1069193>`_.       Science 296 1661-1671.
- Holt et al.             (2002) `The Genome Sequence of the Malaria Mosquito Anophelesd gambiae                                        <http://doi.org/10.1126/science.1076181>`_.       Science 298 129-149.
- Istrail et al.          (2004) `Whole Genome Shotgun Assembly and Comparison of Human Genome Assemblies                               <http://doi.org/10.1073/pnas.0307971100>`_.       PNAS 101 1916-1921.
- Kirkness et al.         (2003) `The Dog Genome: Survey Sequencing and Comparative Analysis                                            <http://doi.org/10.1126/science.1086432>`_.       Science 301 1898-1903.
- Venter et al.           (2004) `Environmental genome shotgun sequencing of the Sargasso Sea                                           <http://doi.org/10.1126/science.1093857>`_.       Science 304 66-74.
- Levy et al.             (2007) `The Diploid Genome Sequence of an Individual Human                                                    <http://doi.org/10.1371/journal.pbio.0050254>`_.  PLoS Biology 0050254
- Rusch et al.            (2007) `The Sorcerer II Global Ocean Sampling Expedition: Northwest Atlantic through Eastern Tropical Pacific <http://doi.org/10.1371/journal.pbio.0050077>`_.  PLoS Biology 1821060.
- Denisov et al.          (2008) `Consensus Generation and Variant Detection by Celera Assembler                                        <http://doi.org/10.1093/bioinformatics/btn074>`_. Bioinformatics 24(8):1035-40
- Miller et al.           (2008) `Aggressive Assembly of Pyrosequencing Reads with Mates                                                <http://doi.org/10.1093/bioinformatics/btn548>`_. Bioinformatics 24(24):2818-2824
- Koren et al.            (2012) `Hybrid error correction and de novo assembly of single-molecule sequencing reads                      <http://doi.org/10.1038/nbt.2280>`_.              Nature Biotechnology, July 2012.
- Koren et al.            (2013) `Reducing assembly complexity of microbial genomes with single-molecule sequencing                     <http://doi.org/10.1186/gb-2013-14-9-r101>`_.     Genome Biology 14:R101.
- Berlin et. al.          (2015) `Assembling Large Genomes with Single-Molecule Sequencing and Locality Sensitive Hashing               <http://doi.org/10.1038/nbt.3238>`_.              Nature Biotechnology. (2015).
- Koren and Walenz et al. (2017) `Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation     <http://doi.org/10.1101/gr.215087.116>`_.         Genome Research. (2017).
