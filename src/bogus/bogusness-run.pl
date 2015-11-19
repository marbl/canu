#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  This file is derived from:
 #
 #    src/AS_BAT/bogusness-run.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2010-DEC-02 to 2013-AUG-01
 #      are Copyright 2010-2011,2013 J. Craig Venter Institute, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz on 2014-DEC-19
 #      are Copyright 2014 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;
use FindBin;

my $prefix = undef;    #  output file name prefix
my $SEQ    = undef;    #  input sequences, fasta
my $REF    = undef;    #  input reference, fasta
my $IDEAL  = undef;    #  input ideal, output from bogus, *.ideal.intervals
my $IGFF3  = undef;    #  input ideal, output from bogus, *.ideal.gff3

#  Assumes this script is NOT installed; it looks to this directory to get
#  jbrowse source and configs.
my $src    = $FindBin::Bin;

#  1.2.1: Needs: (same)
#  1.2:   Needs: Heap::Simple Heap::Simple::Perl Heap::Simple::XS PerlIO::gzip Devel::Size
#  1.1:   Needs: BioPerl, JSON, JSON::XS
#
my $jbrowseVersion = "1.2.1";
my $jbrowse        = "$src/jbrowse-$jbrowseVersion.zip\n";


while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if ($arg eq "-prefix") {
        $prefix = shift @ARGV;

    } elsif ($arg eq "-sequence") {
        $SEQ = shift @ARGV;

    } elsif ($arg eq "-reference") {
        $REF = shift @ARGV;

    } elsif ($arg eq "-ideal") {
        $IDEAL = shift @ARGV;

    } else {
        die "Unknown option '$arg'\n";
    }
}
if (!defined($prefix) || !defined($SEQ) || !defined($REF) || !defined($IDEAL)) {
    print STDERR "usage: $0 -prefix P -sequence A.fasta -reference R.fasta -ideal I.ideal.intervals\n";
    exit(1);
}

if ($IDEAL =~ m/(.*).intervals/) {
    $IDEAL = "$1.intervals";
    $IGFF3 = "$1.gff3"
}
die "IDEAL '$IDEAL' not found.\n"  if (! -e "$IDEAL");
die "IGFF3 '$IGFF3' not found.\n"  if (! -e "$IGFF3");



#
#  Map assembly to reference
#

if (! -e "$prefix.delta") {
    print STDERR "Running nucmer\n";
    my $cmd;

    $cmd .= "nucmer";
    $cmd .= " --maxmatch --coords -p $prefix";
    $cmd .= " $REF";
    $cmd .= " $SEQ";

    system($cmd);
}

if (! -e "$prefix.png") {
    print STDERR "Running mummerplot\n";
    my $cmd;

    $cmd .= "mummerplot";
    $cmd .= " --layout --filter -p $prefix -t png";
    #$cmd .= " --layout          -p $prefix -t png";
    $cmd .= " $prefix.delta";

    system($cmd);
}

#
#  Run bogusness on that mapping
#

if (! -e "$prefix.bogusness.out") {
    print STDERR "Running bogusness\n";
    my $cmd;

    $cmd  = "bogusness \\\n";
    $cmd .= " -reference $REF \\\n";
    $cmd .= " -ideal     $IDEAL \\\n";
    $cmd .= " -nucmer    $prefix.coords \\\n";
    $cmd .= " -output    $prefix \\\n";
    $cmd .= " > $prefix.bogusness.out";

    system($cmd);
}

#
#  Build wiki-friendly output
#

if (! -e "$prefix.bogusness.wiki") {
    my $lastUtg;

    open(F, "sort -k10n -k18n < $prefix.bogusness.out |") or die "Failed to read sort output.\n";
    open(O, "> $prefix.bogusness.wiki") or die "Failed to open '$prefix.bogusness.wiki' for writing.\n";

    print O "{| class=\"wikitable\" border=1\n";
    print O "! unitig !! align num !! utg coords !! gen coords !! status !! ideal type !! ideal index !! ideal coords !! length !! utg cov !! ideal cov !! annotation\n";
    print O "|-\n";

    while (<F>) {
        chomp;

        if (m/^\|\s(utg\d+)\s/) {
            if (!defined($lastUtg)) {
                $lastUtg = $1;
            }
            if ($lastUtg ne $1) {
                print O "| colspan=12 bgcolor=#666666 |\n";
                print O "|-\n";
                $lastUtg = $1;
            }

            s!BEGINSin  \|\| UNIQ!bgcolor="FireBrick" \| BEGINSin  \|\| bgcolor="FireBrick" \| UNIQ!;
            s!ENDSin    \|\| UNIQ!bgcolor="FireBrick" \| ENDSin    \|\| bgcolor="FireBrick" \| UNIQ!;

            s!BEGINSin  \|\| REPT!bgcolor="ForestGreen" \| BEGINSin  \|\| bgcolor="ForestGreen" \| REPT!;
            s!ENDSin    \|\| REPT!bgcolor="ForestGreen" \| ENDSin    \|\| bgcolor="ForestGreen" \| REPT!;

            s!CONTAINED \|\| UNIQ!bgcolor="Indigo" \| CONTAINED \|\| bgcolor="Indigo" \| UNIQ!;

            s!CONTAINED \|\| REPT!bgcolor="SteelBlue" \| CONTAINED \|\| bgcolor="SteelBlue" \| REPT!;

            print O "$_\n";
            print O "|-\n";
        }
    }

    print O "|}\n";

    close(O);
    close(F);
}

#
#  Attempt to create a jbrowse instance.
#
#  cd x010 && tar -xf ../jbrowse.tar && mv jbrowse jbrowse-ctg && cd jbrowse-ctg/ && \
#  ln -s ../bogusness-ctg.gff3 ../bogus.gff3 . && \
#  bin/prepare-refseqs.pl -fasta ../../NC_000913.fasta && \
#  bin/biodb-to-json.pl   --conf bogus.json && \
#  cd ../..
#

system("mkdir $prefix.jbrowse");


open(F, "> $prefix.jbrowse/create.sh") or die "Failed to open '$prefix.jbrowse/create.sh' for writing.\n";
print F "#!/bin/sh\n";
print F "\n";
print F "cd $prefix.jbrowse\n";
print F "\n";
print F "unzip -q $jbrowse\n";
print F "\n";
print F "mv     jbrowse-$jbrowseVersion/* jbrowse/* .\n";
print F "rm -rf jbrowse-$jbrowseVersion\n";
print F "\n";
print F "cp -p $src/bogus-genome.css  genome.css\n";
print F "cp -p $src/bogus-genome.json bogus.json\n";
print F "\n";
print F "ln -s ../$prefix.gff3\n";
print F "\n";
print F "ln -s ../$IGFF3 .\n";
print F "\n";
print F "echo bin/prepare-refseqs.pl -fasta ../$REF\n";
print F "bin/prepare-refseqs.pl -fasta ../$REF\n";
print F "\n";
print F "echo bin/biodb-to-json.pl   --conf bogus.json\n";
print F "bin/biodb-to-json.pl   --conf bogus.json\n";
print F "\n";
close(F);

system("sh $prefix.jbrowse/create.sh");
