#!/usr/bin/perl

use strict;

print STDERR "READING GKPSTORE\n";
my %mate;

open(F, "gatekeeper -dumpfragments -tabular -nouid ../test.gkpStore 2>/dev/null |");
while (<F>) {
    my @v = split '\s+', $_;
    if ($v[3] ne "NA") {
        $mate{$v[1]} = $v[3];
    } else {
        $mate{$v[1]} = 0;
    }
}
close(F);

print STDERR "READING TIGSTORE\n";
my $utg;
my %frgToUnitig;
my %numFragsInUnitig;

open(F, "tigStore -g ../test.gkpStore -t debug.unitigs.aftercontains.tigStore 2 -d layout -U 2>/dev/null |");
while (<F>) {
    if (m/^unitig\s+(\d+)$/) {
        $utg = $1;
    }
    if (m/^FRG.*ident\s+(\d+)\s+/) {
        $frgToUnitig{$1} = $utg;
        $numFragsInUnitig{$utg}++;
    }
}
close(F);


print STDERR "READING UNITIG\n";
$utg = shift @ARGV;

my %utgFrags;

open(F, "tigStore -g ../test.gkpStore -t debug.unitigs.aftercontains.tigStore 2 -d layout -u $utg 2>/dev/null |");
while (<F>) {
    if (m/^FRG.*ident\s+(\d+)\s+.*position\s+(\d+)\s+(\d+)$/) {
        my $fid = $1;
        my $bgn = $2;
        my $end = $3;
        my $mid = $mate{$1};

        $utgFrags{$fid}++;

        next if ($mid == 0);

        my $mtg = $frgToUnitig{$mid};

        next if ($mtg == $utg);

        print STDERR "$1 $bgn,$end -- mate $mid -- in unitig $mtg ($numFragsInUnitig{$mtg} frags)\n";
    }
}
close(F);

