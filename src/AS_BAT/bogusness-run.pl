#!/usr/bin/perl

use strict;

my $prefix = undef;  # "test.004.buildUnitigs";
my $REF    = undef;  #
my $IDEAL  = undef;  # "../../ideal/porphyromonas_gingivalis_w83.flx.fragment.E9T0MN001.ideal.intervals";

$prefix = $ARGV[0]  if (scalar(@ARGV) > 0);
$REF    = $ARGV[1]  if (scalar(@ARGV) > 1);
$IDEAL  = $ARGV[2]  if (scalar(@ARGV) > 2);

die "Must supply 'PREFIX' as first non-option argument.\n"           if (!defined($prefix));
die "Must supply 'REF.fasta' as second non-option argument.\n"       if (!defined($REF));
die "Must supply 'IDEAL.intervals' as third non-option argument.\n"  if (!defined($IDEAL));

if (! -e "$prefix.bogus.out") {
    my $cmd;

    $cmd  = "bogusness \\\n";
    $cmd .= " -reference $REF \\\n";
    $cmd .= " -ideal     $IDEAL \\\n";
    $cmd .= " -nucmer    $prefix.coords \\\n";
    $cmd .= " -output    $prefix \\\n";
    $cmd .= " > $prefix.bogusness.out";

    system($cmd);
}

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
