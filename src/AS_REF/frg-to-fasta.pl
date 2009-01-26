#!/usr/bin/perl
use strict;

while (!eof(STDIN)) {
    $_ = <STDIN>;
    chomp;

    if (m/^acc:(.*)$/) {
        print ">$1\n";
    }

    if (m/^seq:$/) {

        $_ = <STDIN>;
        chomp;

        while ($_ ne ".") {
            print "$_\n";
            $_ = <STDIN>;
            chomp;
        }
    }
}
