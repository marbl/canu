#!/usr/bin/perl

#  Computes the number and cumulative length of regions where two atac
#  mappings disagree.

use strict;

my @alignments;

#  Count the number of small or large differences in location
my @smallH;
my @smallHLen;
my $small = 0;
my $smallLen = 0;
my $large = 0;
my $largeLen = 0;

my $major = 0;
my $majorLen = 0;

my $smallLimit = 40;

open(F, "< overlap.Aannotation");
while (<F>) {
    if (m/^[N!]\s+(\d+):(\d+)-(\d+)\[\s*\d+\].*\s(\d+):\s*(\d+)-\s*(\d+)\).*\s(\d+):\s*(\d+)-\s*(\d+)\)/) {
        my $id1 = $1;
        my $b1  = $2;
        my $e1  = $3;

        my $id2a = $4;
        my $b2a  = $5;
        my $e2a  = $6;

        my $id2b = $7;
        my $b2b  = $8;
        my $e2b  = $9;

        my $diff = 0;

        if ($id2a == $id2b) {
            if ($b2a > $b2b) {
                $diff = $b2a - $b2b;
            } else {
                $diff = $b2b - $b2a;
            }

            if ($diff < $smallLimit) {
                push @alignments, "$1 $2 $3 $4 $5 $6 $7 $8 $9\n";
                $smallH[$diff]++;
                $smallHLen[$diff] += $e1 - $b1;
                $small++;
                $smallLen += $e1 - $b1;
            } else {
                $large++;
                $largeLen += $e1 - $b1;
            }
        } else {
            $major++;
            $majorLen += $e1 - $b1;
        }

        #print STDERR "$_" if ($e1 - $b1 <= 0);
    }
}
close(F);

for (my $i=1; $i<$smallLimit; $i++) {
    if (defined($smallH[$i])) {
        print STDERR "small[ $i ] $smallH[$i]  $smallHLen[$i] bp)\n";
    }
}
print STDERR "small:    $small    ($smallLen bp)\n";
print STDERR "large:    $large    ($largeLen bp)\n";
print STDERR "major:    $major    ($majorLen bp)\n";
