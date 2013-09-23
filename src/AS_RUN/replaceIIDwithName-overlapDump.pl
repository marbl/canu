#!/usr/bin/perl

use strict;

my %IIDtoNAME;

my $fastqUIDmap = shift @ARGV;

print STDERR "Loading UID map from '$fastqUIDmap'.\n";

open (F, "< $fastqUIDmap") or die "Failed to open '$fastqUIDmap'\n";
while (<F>) {
    my @v = split '\s+', $_;

    if      (scalar(@v) == 3) {
        $IIDtoNAME{$v[1]} = $v[2];

    } elsif (scalar(@v) == 6) {
        $IIDtoNAME{$v[1]} = $v[2];
        $IIDtoNAME{$v[4]} = $v[5];

    } else {
        die "unknown format '$_'\n";
    }
}
close(F);


while (<STDIN>) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;

    my @v = split '\s+', $_;

    die "Didn't find IID '$v[0]' in overlap '$_'.\n"  if (!exists($IIDtoNAME{$v[0]}));
    die "Didn't find IID '$v[0]' in overlap '$_'.\n"  if (!exists($IIDtoNAME{$v[1]}));

    $v[0] = $IIDtoNAME{$v[0]};
    $v[1] = $IIDtoNAME{$v[1]};

    print join("\t", @v), "\n";
}
