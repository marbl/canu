#!/usr/bin/perl

use strict;

my %h;

while (<STDIN>) {
    my ($hash, $seq) = split '\s+', $_;
    if (!defined($h{$hash})) {
        $h{$hash} = $seq;
    } else {
        if ($h{$hash} ne $seq) {
            print STDOUT "$hash - $h{$hash}\n";
            print STDOUT "$hash - $seq\n";
        }
    }
}
