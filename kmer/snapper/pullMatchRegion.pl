#!/bin/perl

use lib "/work/assembly/walenzbp/projects/scripts";
use libBri;

while (!eof(STDIN)) {
    my %p = &libBri::readPolish(*STDIN);

    if (defined($p{'raw'})) {
        my $status = "$p{'numMatches'}-$p{'numMatchesN'}-$p{'percentID'}-$p{'matchOrientation'}-$p{'strandPrediction'}";
        my $st = 0;
        my $ed = 0;
        my $ll = 0;

        foreach my $exon (@{$p{'exons'}}) {
            if ($exon->{'GENOMICend'} - $exon->{'GENOMICstart'} > $ll) {
                $st = $p{'dbLo'} + $exon->{'GENOMICstart'};
                $ed = $p{'dbLo'} + $exon->{'GENOMICend'} + 1;
                $ll = $exon->{'GENOMICend'} - $exon->{'GENOMICstart'};
            }
        }

        print "$p{'estID'}\[$p{'estLen'}-$p{'pA'}-$p{'pT'}\] $p{'dbID'}\[$st-$ed\] <$status>\n";
    }
}
