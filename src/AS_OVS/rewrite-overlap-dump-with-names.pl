#!/usr/bin/perl

use strict;

my %IIDtoNAME;

open(F, "< test.gkpStore.fastqUIDmap") or die "Failed to open 'test.gkpStore.fastqUIDmap'\n";
while (<F>) {
    my @v = split '\s+', $_;

    $IIDtoNAME{$v[1]} = $v[2];
}
close(F);


while (<STDIN>) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;

    my @v = split '\s+', $_;

    die if (!exists($IIDtoNAME{$v[0]}));
    die if (!exists($IIDtoNAME{$v[1]}));

    print "$IIDtoNAME{$v[0]}\t$IIDtoNAME{$v[1]}\n";
}
