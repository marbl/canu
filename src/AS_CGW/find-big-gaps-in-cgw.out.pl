#!/usr/bin/perl

#InsertScaffoldContentsIntoScaffold()-- Insert CI    76908   1228bp fwd         0 +- 0                1228 +- 32          was         0 +- 0                1228 +- 32         
#InsertScaffoldContentsIntoScaffold()-- Insert CI   320456  10316bp fwd      2699 +- 53691           13015 +- 53959       was      2699 +- 53691           13015 +- 53959      

my $last = 0;
my $max  = 0;

while (<STDIN>) {
    chomp;

    my @v = split '\s+', $_;

    my $gs = $v[6] - $last;

    $last = $v[6];

    if ($gs > 500000) {
        print "$gs -- $_\n";
    }

    $gapSize[$gs]++;

    $max = ($gs < $max) ? $max : $gs;
}


for (my $i=0; $i<=$max; $i++) {
    if ($gapSize[$i] > 0) {
        print "$i\t$gapSize[$i]\n";
    }
}
