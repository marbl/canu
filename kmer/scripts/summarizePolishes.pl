#!/usr/local/bin/perl

$| = 1;

use strict;

#  We use the libBri package of usefull stuff.  It's located in the same place
#  as this script.  That's what FindBin tells us.
#
use FindBin;
use lib "$FindBin::Bin";
use libBri;

if (scalar @ARGV == 0) {
    print STDERR "usage: $0 <sim4db files to summarize>\n";
    print STDERR "\n";
    print STDERR "NOTE: One summary is generated for all files\n";
    exit;
}

my ($mat, $est, $scf) = &libBri::summarizePolishes(@ARGV);
if ($mat > 0) {
    print "EST-scaffold matches              $mat EST-scaffold pairs ($est different ESTs and $scf scaffolds)\n";
    print "Number of matches per EST         ", $mat / $est, " matches/EST\n";
    print "Number of matches per scaffold    ", $mat / $scf, " matches/scaffold\n";
} else {
    print "EST-scaffold matches              None.\n";
}
print "\n";
