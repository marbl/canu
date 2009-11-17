#!/usr/bin/perl

#  Hack to read two qc files on the command line and report
#  if they differ meaningfully.

use strict;

if (scalar(@ARGV) != 2) {
    die "usage: $0 one.qc two.qc\n";
}

my %gainlossGood;



sub readOrder ($) {
    my $qcname = shift @_;
    my @order;

    open(A, "< $qcname") or die "can't open $qcname\n";
    while (<A>) {
        if (m/^(.*)=(.*)$/) {
            push @order, $1;
        }
    }

    return(@order);
}

sub readQC ($) {
    my $qcname = shift @_;
    my %vals;

    open(A, "< $qcname") or die "can't open $qcname\n";
    while (<A>) {
        if (m/^(.*)=(.*)$/) {
        }

        #  Scaffolds
        if      (m/TotalScaffolds=(\d+)/) {                $vals{'TotalScaffolds'} = $1;
        } elsif (m/TotalContigsInScaffolds=(\d+)/) {       $vals{'TotalContigsInScaffolds'} = $1;
        } elsif (m/TotalBasesInScaffolds=(\d+)/) {         $vals{'TotalBasesInScaffolds'} = $1;

        } elsif (m/N50ScaffoldBases=(\d+)/) {              $vals{'N50ScaffoldBases'} = $1;
        } elsif (m/TotalSpanOfScaffolds=(\d+)/) {          $vals{'TotalSpanOfScaffolds'} = $1;
        }

        if      (m/TotalNumVarRecords=(\d+)/) {            $vals{'TotalNumVarRecords'} = $1;

        #  Contigs
        } elsif (m/N50ContigBases=(\d+)/) {                $vals{'N50ContigBases'} = $1;

        } elsif (m/TotalBigContigs=(\d+)/) {               $vals{'TotalBigContigs'} = $1;
        } elsif (m/BigContigLength=(\d+)/) {               $vals{'BigContigLength'} = $1;
        } elsif (m/BigContigsPercentBases=(\d+.\d+)/) {    $vals{'BigContigsPercentBases'} = $1;

        } elsif (m/TotalSmallContigs=(\d+)/) {             $vals{'TotalSmallContigs'} = $1;
        } elsif (m/SmallContigLength=(\d+)/) {             $vals{'SmallContigLength'} = $1;
        } elsif (m/SmallContigsPercentBases=(\d+.\d+)/) {  $vals{'SmallContigsPercentBases'} = $1;

        } elsif (m/TotalDegenContigs=(\d+)/) {             $vals{'TotalDegenContigs'} = $1;
        } elsif (m/DegenContigLength=(\d+)/) {             $vals{'DegenContigLength'} = $1;
        } elsif (m/DegenPercentBases=(\d+.\d+)/) {         $vals{'DegenPercentBases'} = $1;

        } elsif (m/TotalUUnitigs=(\d+)/) {                 $vals{'TotalUUnitigs'} = $1;
        }

        #  Surrogates
        if      (m/TotalSurrogates=(\d+)/) {               $vals{'TotalSurrogates'} = $1;
        } elsif (m/SurrogateInstances=(\d+)/) {            $vals{'SurrogateInstances'} = $1;
        } elsif (m/SurrogateLength=(\d+)/) {               $vals{'SurrogateLength'} = $1;
        } elsif (m/SurrogateInstanceLength=(\d+)/) {       $vals{'SurrogateInstanceLength'} = $1;
        }

        #  Reads
        if      (m/ReadsWithNoMate=(\d+)/) {               $vals{'ReadsWithNoMate'} = $1;
        } elsif (m/ReadsWithGoodMate=(\d+)/) {             $vals{'ReadsWithGoodMate'} = $1;
        } elsif (m/ReadsWithBadShortMate=(\d+)/) {         $vals{'ReadsWithBadShortMate'} = $1;
        } elsif (m/ReadsWithBadLongMate=(\d+)/) {          $vals{'ReadsWithBadLongMate'} = $1;

        } elsif (m/TotalUsableReads=(\d+)/) {              $vals{'TotalUsableReads'} = $1;
        } elsif (m/AvgClearRange=(\d+.\d+)/) {             $vals{'AvgClearRange'} = $1;
        } elsif (m/ContigReads=(\d+)/) {                   $vals{'ContigReads'} = $1;
        } elsif (m/BigContigReads=(\d+)/) {                $vals{'BigContigReads'} = $1;
        } elsif (m/SmallContigReads=(\d+)/) {              $vals{'SmallContigReads'} = $1;
        } elsif (m/DegenContigReads=(\d+)/) {              $vals{'DegenContigReads'} = $1;
        } elsif (m/PlacedSurrogateReads=(\d+)/) {          $vals{'PlacedSurrogateReads'} = $1;
        } elsif (m/SingletonReads=(\d+)/) {                $vals{'SingletonReads'} = $1;
        } elsif (m/ChaffReads=(\d+)/) {                    $vals{'ChaffReads'} = $1;
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

    foreach my $k (readOrder($ARGV[0])) {
        my $d =  $b{$k} - $a{$k};
        my $m = ($b{$k} > $a{$k}) ? $b{$k} : $a{$k};

        $a{$k} = int($a{$k} * 1000) / 1000;
        $b{$k} = int($b{$k} * 1000) / 1000;
        $d     = int($d     * 1000) / 1000;

        if ($m > 0) {
            my $p = int(1000 * $d / $m) / 10;

            my $K = substr("$k                    ", 0, 20);
            my $A = substr("$a{$k}                    ", 0, 20);
            my $B = substr("$b{$k}                    ", 0, 20);

            if (abs($d / $m) > 0.001) {
                $p = "+$p" if ($p > 0.0);
                print STDOUT "$K\t$p%\t$d\tA:$A\tB:$B\n";
            }
        }
    }
}

