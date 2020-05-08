#!/usr/bin/env perl

use strict;

my $tigid;
my $readlist;

while (<STDIN>) {
    if (m/^tig\s(\d+)/) {
        $tigid = $1;
    }

    if (m/^read\s+(\d+)\s+anchor.*position\s+(\d+)\s+(\d+)\s/) {
        my $bgn = $2;
        my $end = $3;
        my $rid = $1 . (($bgn < $end) ? "+" : "-");

        if (defined($readlist)) {
            $readlist .= ",$rid";
        } else {
            $readlist  = $rid;
        }
    }

    if (m/^tigend/) {
        print "tig$tigid $readlist\n";

        undef $tigid;
        undef $readlist;
    }
}
