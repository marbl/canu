#!/usr/bin/perl

#  Reads a list of numbers on stdin, computes the n50.
#
#  If there are two numbers per line, they are assumed to be
#  a begin-end pair.

#  grep "M u " ATAC/atac.shift.atac | cut -d' ' -f 7 | perl n50.pl 3076782067
#  grep "M u " ATAC/box2.shift.atac | cut -d' ' -f 7 | perl n50.pl 3076782067

my @values;

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

    push @values, $val;
}

if (scalar(@ARGV) > 0) {
    $totalLength = int($ARGV[0]);
} else {
    foreach my $v (@values) {
        $totalLength += $v;
    }
}

@values = sort { $a <=> $b } @values;

for (my $nvalue = 1; $nvalue <= 100; $nvalue += 1) {
    my $limit = $nvalue * $totalLength / 100;
    my $iter  = 0;
    my $sum   = 0;

    while (($sum < $limit) && ($iter < scalar(@values))) {
        $sum += $values[$iter++];
    }

    print STDOUT "$nvalue $limit : $values[$iter-1]\n";
}
