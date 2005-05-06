#!/usr/bin/perl

#  Reads a list of numbers on stdin, computes a (blocked) histogram.
#
#  If there are two numbers per line, they are assumed to be
#  a begin-end pair.

#  grep "M u " ATAC/atac.shift.atac | cut -d' ' -f 7 | perl run-length-histogram.pl > atac.histogram
#  grep "M u " ATAC/box2.shift.atac | cut -d' ' -f 7 | perl run-length-histogram.pl > box2.histogram

my @histogram;
my $blocksize = 1000;

while (<STDIN>) {
    s/^\s+//;
    s/\s$//;
    s/\s+/ /;

    my @vals = split '\s+', $_;
    my $val;

    if      (scalar(@vals) == 1) {
        $val = $vals[0];
    } elsif (scalar(@vals) == 1) {
        $val = $vals[1] - $vals[0];
        $val = $vals[0] - $vals[1] if ($val < 0);
    } else {
    }

    $val = ($val / $blocksize);
    $histogram[$val]++;
}

my $max = scalar(@histogram) + 1;
for (my $i=0; $i<$max; $i++) {
    $histogram[$i] = 0 if ($histogram[$i] == 0);
    print "$i $histogram[$i]\n";
}

