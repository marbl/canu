#!/usr/bin/perl

use strict;

my $usr = 0;
my $sys = 0;
my $wal = 0;

foreach my $f (@ARGV) {
    open(F, "< $f");
    while (<F>) {
        if (m/userTime:\s+(\d+.\d)/)   { $usr += $1; }
        if (m/systemTime:\s+(\d+.\d)/) { $sys += $1; }
        if (m/wallTime:\s+(\d+.\d)/)   { $wal += $1; }
    }
    close(F);
}

print "user: $usr   system: $sys   wall: $wal\n";
