#!/usr/local/bin/perl

#
#
#

my $currentestid = -1;
my $outprefix;
my @outline;
my $num          = 0;

my $maxscore = 0;

while ($currentestid < 222439) {
    $_ = <STDIN>;

    my ($dir, $junk, $estid, $junk, $chr, $beg, $end, $junk, $s1, $s2, $sm, $junk, $i, $c) = split '\s+', $_;

    if ($currentestid == $estid) {
        push @outline, "$s1,$s2,$i,$c";
        $maxscore = $sm;
        $num++;
    } else {
        @outline = sort { $b <=> $a } @outline;

        if (defined($outprefix)) {
            print "$outprefix\t$maxscore\t$num\t";

            my $a = $,;
            my $b = $\;
            $, = " ";
            $\ = "\n";

            print @outline;

            $, = $a;
            $\ = $b;
        }

        $currentestid = $estid;
        $outprefix    = "$estid";

        undef @outline;

        push @outline, "$s1,$s2,$i,$c";
        $maxscore = $sm;
        $num      = 1;
    }
}

