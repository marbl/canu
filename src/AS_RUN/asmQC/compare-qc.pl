#!/usr/bin/perl

#  Hack to read two qc files on the command line and report
#  if they differ meaningfully.

use strict;

if (scalar(@ARGV) != 2) {
    die "usage: $0 one.qc two.qc\n";
}

sub readQC ($) {
    my $qcname = shift @_;
    my %vals;

    open(A, "< $qcname") or die "can't open $qcname\n";
    while (<A>) {
        if      (m/TotalScaffolds=(\d+)/) {
            $vals{'TotalScaffolds'} = $1;
        } elsif (m/TotalContigsInScaffolds=(\d+)/) {
            $vals{'TotalContigsInScaffolds'} = $1;
        } elsif (m/TotalBasesInScaffolds=(\d+)/) {
            $vals{'TotalBasesInScaffolds'} = $1;
        } elsif (m/N50ScaffoldBases=(\d+)/) {
            $vals{'N50ScaffoldBases'} = $1;
        } elsif (m/TotalSpanOfScaffolds=(\d+)/) {
            $vals{'TotalSpanOfScaffolds'} = $1;
        } elsif (m/TotalNumVarRecords=(\d+)/) {
            $vals{'TotalNumVarRecords'} = $1;
        } elsif (m/N50ContigBases=(\d+)/) {
            $vals{'N50ContigBases'} = $1;
        } elsif (m/TotalBigContigs=(\d+)/) {
            $vals{'TotalBigContigs'} = $1;
        } elsif (m/BigContigLength=(\d+)/) {
            $vals{'BigContigLength'} = $1;
        } elsif (m/TotalSmallContigs=(\d+)/) {
            $vals{'TotalSmallContigs'} = $1;
        } elsif (m/SmallContigLength=(\d+)/) {
            $vals{'SmallContigLength'} = $1;
        } elsif (m/SmallContigsPercentBases=(\d+.\d+)/) {
            $vals{'SmallContigsPercentBases'} = $1;
        } elsif (m/TotalDegenContigs=(\d+)/) {
            $vals{'TotalDegenContigs'} = $1;
        } elsif (m/DegenContigLength=(\d+)/) {
            $vals{'DegenContigLength'} = $1;
        } elsif (m/DegenPercentBases=(\d+.\d+)/) {
            $vals{'DegenPercentBases'} = $1;
        } elsif (m/TotalUUnitigs=(\d+)/) {
            $vals{'TotalUUnitigs'} = $1;
        } elsif (m/TotalSurrogates=(\d+)/) {
            $vals{'TotalSurrogates'} = $1;
        } elsif (m/SurrogateInstances=(\d+)/) {
            $vals{'SurrogateInstances'} = $1;
        } elsif (m/SurrogateLength=(\d+)/) {
            $vals{'SurrogateLength'} = $1;
        } elsif (m/SurrogateInstanceLength=(\d+)/) {
            $vals{'SurrogateInstanceLength'} = $1;
        } elsif (m/ReadsWithNoMate=(\d+)/) {
            $vals{'ReadsWithNoMate'} = $1;
        } elsif (m/ReadsWithGoodMate=(\d+)/) {
            $vals{'ReadsWithGoodMate'} = $1;
        } elsif (m/ReadsWithBadShortMate=(\d+)/) {
            $vals{'ReadsWithBadShortMate'} = $1;
        } elsif (m/ReadsWithBadLongMate=(\d+)/) {
            $vals{'ReadsWithBadLongMate'} = $1;
        } elsif (m/TotalReads=(\d+)/) {
            $vals{'TotalReads'} = $1;
        } elsif (m/AvgClearRange=(\d+.\d+)/) {
            $vals{'AvgClearRange'} = $1;
        } elsif (m/ReadsInContigs=(\d+)/) {
            $vals{'ReadsInContigs'} = $1;
        } elsif (m/PlacedSurrogateReads=(\d+)/) {
            $vals{'PlacedSurrogateReads'} = $1;
        } elsif (m/ContigsOnly=(\d+)/) {
            $vals{'ContigsOnly'} = $1;
        } elsif (m/Content=(\d+.\d+)/) {
            $vals{'Content'} = $1;
        }
    }
    close(A);

    return(%vals);
}

my %a = readQC($ARGV[0]);
my %b = readQC($ARGV[1]);

my $doOutput = 0;

foreach my $k (keys %a) {
    my $d = $a{$k} - $b{$k};
    my $m = ($a{$k} > $b{$k}) ? $a{$k} : $b{$k};
    if ($m > 0) {
        my $p = int(1000 * $d / $m) / 10;
        my $K = substr("$k            ", 0, 20);
        if (abs($d / $m) > 0.01) {
            $doOutput++;
        }
    }
}

if ($doOutput) {
    print STDOUT "----------------------------------------\n";
    print STDOUT "A:$ARGV[0]\n";
    print STDOUT "B:$ARGV[1]\n";

    foreach my $k (keys %a) {
        my $d = $a{$k} - $b{$k};
        my $m = ($a{$k} > $b{$k}) ? $a{$k} : $b{$k};
        if ($m > 0) {
            my $p = int(1000 * $d / $m) / 10;
            my $K = substr("$k            ", 0, 20);
            if (abs($d / $m) > 0.01) {
                $p = "+$p" if ($p > 0.0);
                print STDOUT "$K\t$p%\tA:$a{$k}\tB:$b{$k}\n";
            }
        }
    }
}

