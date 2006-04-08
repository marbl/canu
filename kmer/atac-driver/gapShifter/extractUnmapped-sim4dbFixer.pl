#!/usr/bin/perl

#  Fixes up the IID's and coords of the sim4db mapped regions from extractUnmapped.
#  Reads polishes from stdin, writes fixed polishes to stdout.

use strict;
use lib "/bioinfo/assembly/walenz/src/genomics/scripts";
use sim4polish;
$| = 1;

while (!eof(STDIN)) {
    my %p = &sim4polish::readPolish(*STDIN);

    if ($p{'raw'}) {
        my $estIID;
        my $estBeg = 0;
        my $estEnd = 0;

        my $dbIID;
        my $dbBeg = 0;
        my $dbEnd = 0;

        if ($p{'estDefLine'} =~ m/extracted\s+from\s+iid\s+(\d+)\s+pos\s+(\d+)\s+(\d+)\s+/) {
            $estIID = $1;
            $estBeg = $2;
            $estEnd = $3;
        }

        if ($p{'dbDefLine'} =~ m/extracted\s+from\s+iid\s+(\d+)\s+pos\s+(\d+)\s+(\d+)\s+/) {
            $dbIID = $1;
            $dbBeg = $2;
            $dbEnd = $3;
        }

        if (defined($estIID)) {
            $p{'estID'} = $estIID;
            foreach my $exon (@{@p{'exons'}}) {
                $exon->{'cDNAstart'} += $estBeg;
                $exon->{'cDNAend'}   += $estBeg;
            }
        }
        if (defined($dbIID)) {
            $p{'dbID'} = $dbIID;
            foreach my $exon (@{@p{'exons'}}) {
                $exon->{'GENOMICstart'} += $dbBeg;
                $exon->{'GENOMICend'}   += $dbBeg;
            }
        }

        #  normalize
        foreach my $exon (@{@p{'exons'}}) {
            $exon->{'GENOMICstart'} += $p{'dbLo'};
            $exon->{'GENOMICend'}   += $p{'dbLo'};
        }

        $p{'dbLo'}   = 0;
        $p{'dbHi'}   = 0;
        $p{'estLen'} = 0;

        $p{'raw'} = &sim4polish::updatePolish(%p);
        print $p{'raw'};
    }
}
