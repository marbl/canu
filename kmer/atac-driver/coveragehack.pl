#!/usr/bin/perl

#  Examines an ATAC file for overlaps in the first assembly.

use strict;
my $lastBeg = 0;
my $lastEnd = 0;
my $thisBeg = 0;
my $thisEnd = 0;

open(F, "grep B34LC:0 boxfiller_out_v2.blocks_covLT.8_filtered.matches.atac | sort -k6n |");
while (<F>) {
    print $_;
    my @vals = split '\s+', $_;

    $thisBeg = $vals[5];
    $thisEnd = $vals[5] + $vals[6];

    if ($lastEnd > $thisBeg) {
        printf("Overlap by %d -- last is [%d-%d] this is [%d-%d]\n",
               $lastEnd - $thisBeg,
               $lastBeg, $lastEnd, $thisBeg, $thisEnd);
    }

    $lastBeg = $thisBeg;
    $lastEnd = $thisEnd;
}
close(F);
