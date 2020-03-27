#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

use strict;

#  Generate a circular 'genome' with nUnique pieces and nUnique+1 repeats, the pair at the start/end
#  forming the circularization.  Assemble it, and compare against 'reference'.

my $repeatSize =   10000;
my $uniqueSize =   90000;
my $nUnique    =       4;

my $readLenL   = 4000;  #  Long reads
my $readLenM   = 3500;  #  Medium reads
my $readLenS   = 3000;  #  Short reads

my $coverageL  = 30;
my $coverageM  = 1;
my $coverageS  = 1;

sub randomSequence ($) {
    my $L = shift @_;

    open(F, "leaff -G 1 $L $L |");
    $_ = <F>;  chomp;
    $_ = <F>;  chomp;
    close(F);

    return($_);
}

if (! -e "reads.genome.fasta") {
    my $len = 0;

    my $R    = randomSequence($repeatSize);
    my $Rlen = length($R);

    open(O, "> reads.genome.fasta");

    print O ">G\n";

    for (my $r=0; $r < $nUnique; $r++) {
        print STDERR "REPEAT $len ", $len + $Rlen, "\n";
        print O "$R\n";
        $len += $Rlen;

        my $U    = randomSequence($uniqueSize);
        my $Ulen = length($U);

        print STDERR "UNIQUE $len ", $len + $Ulen, "\n";
        print O "$U\n";
        $len += $Ulen;
    }

    print STDERR "REPEAT $len ", $len + $Rlen, "\n";
    print O "$R\n";
    $len += $Rlen;

    close(O);
}


if (! -e "reads.s.fastq") {
    system("fastqSimulate -f reads.genome.fasta -o readsL -l $readLenL -x $coverageL -em 0.005 -ei 0 -ed 0 -se");
    system("fastqSimulate -f reads.genome.fasta -o readsM -l $readLenM -x $coverageM -em 0.005 -ei 0 -ed 0 -se");
    system("fastqSimulate -f reads.genome.fasta -o readsS -l $readLenS -x $coverageS -em 0.005 -ei 0 -ed 0 -se");
}

my $canu;

$canu  = "canu -assemble -p test -d test";
$canu .= " ovlMerThreshold=0";
$canu .= " useGrid=0 genomeSize=1m ovlThreads=24";
$canu .= " -pacbio-corrected readsL.s.fastq";
$canu .= " -pacbio-corrected readsM.s.fastq";
$canu .= " -pacbio-corrected readsS.s.fastq";

system($canu);

system("dotplot.sh UTG reads.genome.fasta test/test.unitigs.fasta");
system("dotplot.sh CTG reads.genome.fasta test/test.contigs.fasta");

exit(0);

