#!/usr/local/bin/perl

$| = 1;
use strict;
use lib "/work/assembly/walenzbp/projects/scripts";
use libBri;


#  Takes a fasta file and a set of polishes, and rewrites them both,
#  removing sequences that aren't used.
#
#  For example, if the polishes only use cDNA #392 out of 4000, then
#  the resulting fasta will contain just cDNA #392, and the index in
#  the output polishes file will be correct.
#

if (scalar(@ARGV) != 4) {
    print STDERR "usage: $0 in.polishes in.fasta out.polishes out.fasta\n";
    exit;
}

die "Can't find $ARGV[0]\n" if (! -e $ARGV[0]);
die "Can't find $ARGV[1]\n" if (! -e $ARGV[1]);
#die "File exists: $ARGV[2]\n" if (-e $ARGV[2]);
#die "File exists: $ARGV[3]\n" if (-e $ARGV[3]);

#
#  Read the polishes
#

my %IDXs;
my $idx = 0;

open(F, "< $ARGV[0]");
open(O, "> $ARGV[2]");
while (!eof(F)) {
    my %p = &libBri::readPolish(*F);
    if (defined($p{'raw'})) {
        if (!defined($IDXs{$p{'estID'}})) {
            $IDXs{$p{'estID'}} = $idx;
            $idx++;
        }

        $p{'estID'} = $IDXs{$p{'estID'}};
        print O &libBri::updatePolishInfoLine(%p);
    }
}
close(F);
close(O);

print STDERR "Calling leaff for " . scalar(keys %IDXs) . " fasta sequences.\n";


#
#  Use leaff to rewrite the fasta
#
open(O, "> junk.$$.reducePolishes");
foreach my $k (sort {$a <=> $b} keys %IDXs) {
    print O "$k\n";
}
close(O);

system("/work/assembly/walenzbp/projects/releases/leaff -F $ARGV[1] -q junk.$$.reducePolishes > $ARGV[3]");

unlink "junk.$$.reducePolishes";
