#!/usr/bin/perl

use strict;

my @IIDtoNAME;

open(F, "< test.gkpStore.fastqUIDmap") or die "failed to open test.gkpStore.fastqUIDmap\n";
while (<F>) {
    my ($uid, $iid, $name) = split '\s+', $_;

    $IIDtoNAME[$iid] = $name;
}
close(F);

open(F, "overlapStore -d5 -d3 -d test.obtStore |");
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    print "$IIDtoNAME[$v[0]]\t$IIDtoNAME[$v[1]]\t$v[2]\t$v[3]\t$v[4]\t$v[5]\t$v[6]\t$v[7]\n";
}
close(F);
