
This directory contains a variety of utilities that I found to be necessary when
performing the initial analyses of the 454 data.

The central utilities used for shredding 454 assemblies into 3730 pseudoreads are
the following:

	1. Generate_NonShallow_Contigs.pl
	2. Shred_Contigs.pl
	3. FASTA_to_frg_file.pl

These will respectively, extract contigs from an ace file, shred them, then convert them
into a frag file for the Celera assembler.

The other utilities were for correcting unanticipated sequence data in the 454 output such as
finding N's in the reads, or duplicate reads in the same contig.  There are also utilities
for generating/extracting statistics, such as depth, from the defline of the 454 fasta files.

- Kelvin Li (January 2006)

