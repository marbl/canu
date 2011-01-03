#!/usr/bin/perl

use strict;

my $acc = "";

while (!eof(STDIN)) {
    $_ = <STDIN>;
    chomp;

    if (m/^acc:(.*)$/) {
        $acc = $1;
    }

    if (m/^seq:$/) {
        print ">$acc\n";

        $_ = <STDIN>;
        chomp;

        while ($_ ne ".") {
            print "$_\n";
            $_ = <STDIN>;
            chomp;
        }
    }
}
