#!/usr/bin/perl

use strict;

my @logDepth;

for (my $i=5; $i<20; $i++) {
    $logDepth[$i] = "$i";
}

my $inName = "unitigger.err";
if (scalar(@ARGV) > 0) {
    $inName = shift @ARGV;
}

open(F, "< $inName");
while (<F>) {
    if (m/checkUnitigMembership\(\)--\s+(\d+)\s+\(\s+\d+-\s+\d+\)\s+(\d+)/) {
        $logDepth[$1] .= "\t$2";
    }
}
close(F);

for (my $i=5; $i<20; $i++) {
    print "$logDepth[$i]\n";
}
