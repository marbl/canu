#!/usr/bin/perl

use strict;

#  Computes the number and cumulative length of regions where two atac
#  mappings disagree.  
#
#  Reports when the region maps to:
#    small -- same sequence, close together
#    large -- same scaffold, not close together
#    major -- different scaffold
#
#  Automagically generates a plot
#

if (scalar(@ARGV != 2)) {
    print STDERR "usage: $0 some.atac outprefix\n";
    exit(1);
}

my $filename   = shift @ARGV;
my $outprefix  = shift @ARGV;
my $smallLimit = 400;

my @smallH;
my @smallHLen;
my $large = 0;
my $largeLen = 0;
my $major = 0;
my $majorLen = 0;

open(F, "< $filename") or die "Failed to open $filename.\n";
while (<F>) {
    if (m/^[N!]\s+(\d+):(\d+)-(\d+)\[\s*\d+\].*\s(\d+):\s*(\d+)-\s*(\d+)\).*\s(\d+):\s*(\d+)-\s*(\d+)\)/) {
        my ($id1,  $b1,  $e1)  = ($1, $2, $3);
        my ($id2a, $b2a, $e2a) = ($4, $5, $6);
        my ($id2b, $b2b, $e2b) = ($7, $8, $9);

        if ($id2a == $id2b) {
            my $diff;
            $diff = $b2b - $b2a;
            $diff = $b2a - $b2b if ($b2a > $b2b);

            if ($diff < $smallLimit) {
                $smallH[$diff]++;
                $smallHLen[$diff] += $e1 - $b1;
            } else {
                $large++;
                $largeLen += $e1 - $b1;
            }
        } else {
            $major++;
            $majorLen += $e1 - $b1;
        }
    }
}
close(F);

#  output is
#
#  distance away
#  number of regions
#  number of bp in those regions
#  cumulative number of regions
#  cumulative number of bp in those regions
#


my $sumH    = 0;
my $sumHLen = 0;
open(F, "> $outprefix.dat");
for (my $i=1; $i<$smallLimit; $i++) {
    $sumH    += $smallH[$i];
    $sumHLen += $smallHLen[$i] / 10;
    print F "$i $smallH[$i] $smallHLen[$i] $sumH $sumHLen\n" if (defined($smallH[$i]));
}
close(F);

print STDERR "at most  $smallLimit bp away: $sumH regions $sumHLen bp\n";
print STDERR "at least $smallLimit bp away: $large regions $largeLen bp\n";
print STDERR "different sequence: $major regions $majorLen bp\n";

open(F, "> $outprefix.gnuplot");
print F "set terminal postscript color\n";
print F "set output \"$outprefix.ps\"\n";
print F "set xlabel \"bp Difference in Match Location\"\n";
print F "set ylabel \"\"\n";
print F "plot [][0:300000] \"$outprefix.dat\" using 2 with lines title \"Number of Regions\", \\\n";
print F "                  \"$outprefix.dat\" using 3 with lines title \"bp in Regions\", \\\n";
print F "                  \"$outprefix.dat\" using 4 with lines title \"Cumulative Number of Regions\", \\\n";
print F "                  \"$outprefix.dat\" using 5 with lines title \"Cumulative bp in Regions / 10\"\n";
print F "plot [0:100][0:300000] \"$outprefix.dat\" using 2 with lines title \"Number of Regions\", \\\n";
print F "                  \"$outprefix.dat\" using 3 with lines title \"bp in Regions\", \\\n";
print F "                  \"$outprefix.dat\" using 4 with lines title \"Cumulative Number of Regions\", \\\n";
print F "                  \"$outprefix.dat\" using 5 with lines title \"Cumulative bp in Regions / 10\"\n";
close(F);

system("gnuplot < $outprefix.gnuplot");


