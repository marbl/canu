#!/usr/bin/perl

use strict;

my $prefix = undef;    #  output file name prefix
my $SEQ    = undef;    #  input sequences, fasta
my $REF    = undef;    #  input reference, fasta
my $IDEAL  = undef;    #  input ideal, output from bogus, *.ideal.intervals

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

#
#  Map assembly to reference
#

print STDERR "Running nucmer\n";
if (! -e "$prefix.delta") {
    my $cmd;

    $cmd .= "nucmer";
    $cmd .= " --maxmatch --coords -p $prefix";
    $cmd .= " $REF";
    $cmd .= " $SEQ";

    system($cmd);
}

if (! -e "$prefix.png") {
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

exit;

#
#  Attempt to create a jbrowse instance.
#
#  cd x010 && tar -xf ../jbrowse.tar && mv jbrowse jbrowse-ctg && cd jbrowse-ctg/ && \
#  ln -s ../bogusness-ctg.gff3 ../bogus.gff3 . && \
#  bin/prepare-refseqs.pl -fasta ../../NC_000913.fasta && \
#  bin/biodb-to-json.pl   --conf bogus.json && \
#  cd ../..
#

#open(F, "> $prefix.jbrowse.sh") or die "Failed to open '$prefix.jbrowse.sh' for writing.\n";
#print F "#!/bin/sh\n";
#print F "cd $prefix\n";
#print F "tar -xf ../jbrowse.tar\n";
#print F "mv jbrowse/* .\n";
#print F "rmdir jbrowse\n";
#print f "ln -s ../bogusness-ctg.gff3 
