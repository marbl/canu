#!/usr/bin/perl

use strict;

my $line;

open(F, "< the-plan");
while (<F>) {
    chomp;
    s/\s+$//;

    if (m/^\S/) {
        if (!defined($line)) {
            $line = $_;
        } else {
            $line .= " $_";
        }
    } else {
        print "$line\n";
        print "$_\n";
        undef $line;
    }
}
close(F);

