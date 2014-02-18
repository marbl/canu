#!/usr/bin/perl

use strict;

#  Checks cgw.out for large gaps.  Only checks the two scaffolds that are merged.  Does not check the newly formed scaffold.

my $lasbgn = 0;
my $lasend = 0;
my $lasctg = 0;
my $dist   = 0;
my $first  = undef;

my $scf    = 0;

my $lastline = undef;
my $thisline = undef;

while (<STDIN>) {
    $lastline = $thisline;
    $thisline = $_;

    if (m/InsertScaffoldContentsIntoScaffold\(\)--\s+Insert\s+scaffold\s(\d+)/) {
        $lasbgn = 0;
        $lasend = 0;
        $lasctg = 0;
        $dist   = 0;
        $first  = undef;
        $scf    = $1;
    }

    if (m/InsertScaffoldContentsIntoScaffold\(\)--\s+Insert\s+CI\s+(\d+)\s+(\d+bp)\s+(.*)\s+(\d+)\s+..\s+(\d+)\s+(\d+)\s+..\s+(\d+)\s+was\s+(\d+)\s+..\s+(\d+)\s+(\d+)\s+..\s+(\d+)\s+$/) {
        my $thsctg = $1;
        my $thsbgn = ($4 < $6) ? $4 : $6;
        my $thsend = ($4 < $6) ? $6 : $4;

        #                             coords increasing     coords decreasing
        $dist = ($lasbgn < $thsbgn) ? ($thsbgn - $lasend) : ($lasbgn - $thsend);

        if (defined($first)) {
            if ($dist > 150000) {
                print "GAP:\t$dist\tCI\t$lasctg\tpos\t$lasbgn\t$lasend\tto\tCI\t$thsctg\tpos\t$thsbgn\t$thsend\tscf\t$scf\n";
                print $lastline;
                print $thisline;

            } else {
                print "GAP:\t$dist\tCI\t$lasctg\tpos\t$lasbgn\t$lasend\tto\tCI\t$thsctg\tpos\t$thsbgn\t$thsend\tscf\t$scf\n";
            }
        }

        $lasbgn  = $thsbgn;
        $lasend  = $thsend;
        $lasctg  = $thsctg;
        $first  .= $_;
    } else {
    }
}
