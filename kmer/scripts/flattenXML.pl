#!/usr/local/bin/perl

use FileHandle;
use strict;

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

#
#  Flattens a heirarchy of XML into one single XML file.
#  It will find all the files in the directories.  Output is on stdout.
#

if (scalar @ARGV < 1) {
    print STDOUT "usage: $0 outprefix indir indir indir ...\n";
    exit;
}

my $outPrefix = shift @ARGV;
my @inp;

#
#  Find all the input files, and store them sorted by scaffold.
#
print STDERR "Finding XML files.\n";
foreach my $i (@ARGV) {
    open(F, "find $i -type f -name '*.gbf' -print |");
    while (<F>) {
        chomp;
        push @inp, $_;
    }
    close(F);
}

#
#  Check that all the inputs are the same assembly
#
print STDERR "Finding assembly version / taxon\n";
my %check;
my $t;
my $av;
my %outF;
my $fh;
foreach my $f (@inp) {
    open(F, "< $f");
    $t = <F>;
    chomp $t;

    if ($t =~ m/assembly_version="(\d+)"/) {
        $av = $1;
    } else {
        print STDERR "Can't find assembly_version in:\n";
        print STDERR "  '$t'\n";
        exit;
    }

    if (!defined($check{$t})) {
        print STDERR "$outPrefix-$av.gbf\n";

        $outF{$t} = new FileHandle;
        $outF{$t}->open("> $outPrefix-$av.gbf");

        $fh = $outF{$t};
        print $fh "$t\n";

        $check{$t} = 1;
    }
}

#
#  For each scaffold we found, write new output
#
print STDERR "Flattening\n";
foreach my $f (@inp) {
    #print STDERR "$f\r";

    open(F, "< $f");
    $_  = <F>;
    chomp;
    $fh = $outF{$_};
    while (<F>) {
        if ($_ ne "</game>\n") {
            print $fh $_;
        }
    }
    close(F);
}

foreach my $t (keys %check) {
    $fh = $outF{$t};
    print $fh "</game>\n";
    $outF{$t}->close;
}
