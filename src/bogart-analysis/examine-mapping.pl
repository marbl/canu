#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

my $TOOSHORT = 500;

my $FILE     = "test.004.buildUnitigs";
$FILE = shift @ARGV if (scalar(@ARGV) > 0);

if ((! -e "$FILE.tigStore") ||
    (! -e "$FILE.fasta")) {
    die "Missing tigStore or fasta.  Run build-fasta.pl\n";
}

########################################

if (! -e "$FILE.coords") {
    my $cmd;

    $cmd  = "nucmer";
    $cmd .= " --maxmatch --coords -p $FILE";
    $cmd .= " /work/FRAGS/porphyromonas_gingivalis_w83/reference/AE015924.fasta";
    $cmd .= " $FILE.fasta";

    system($cmd);
}

########################################
#  Load the lengh of each sequence.

my %length;

if (-e "$FILE.fasta") {
    open(F, "< $FILE.fasta") or die;
    while (<F>) {
        if (m/^>(utg\d+)\s+len=(\d+)$/) {
            $length{$1} = $2;
        }
    }
    close(F);
} else {
    print STDERR "No $FILE.fasta file found, cannot .....\n";
}

########################################
#  Pass one, load the coords.  Analyze anything with one match, counting
#  the amount of coverage.

my %nucmer;

if (-e "$FILE.coords") {
    open(F, "< $FILE.coords") or die;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    while (<F>) {
        s/^\s+//;
        s/\s$+//;
        my @v = split '\s+', $_;

        my $utgBgn = ($v[3] < $v[4]) ? $v[3] : $v[4];
        my $utgEnd = ($v[3] < $v[4]) ? $v[4] : $v[3];
        my $genBgn = ($v[0] < $v[1]) ? $v[0] : $v[1];
        my $genEnd = ($v[0] < $v[1]) ? $v[1] : $v[0];

        #  Rearrange the coords so that bgn is ALWAYS less than end (we lose orientation).
        my $str = "$utgBgn\t$utgEnd\t$genBgn\t$genEnd";

        if (exists($nucmer{$v[12]})) {
            $nucmer{$v[12]} .= "\n" . $str;
        } else {
            $nucmer{$v[12]} = $str;
        }
    }
    close(F);
} else {
    die "No $FILE.coords??\n";
}

########################################
#  For things with one match, report spurs.

my $spur = 0;
my $spurShort = 0;
my $confirmed = 0;
my $confirmedBP = 0;

foreach my $utg (keys %nucmer) {
    my @m = split '\n', $nucmer{$utg};

    if (scalar(@m) > 1) {
        @m = sort {$a <=> $b} @m;
        $nucmer{$utg} = join "\n", @m;
        next;
    }

    my ($utgBgn, $utgEnd, $genBgn, $genEnd) = split '\s+', $nucmer{$utg};

    my $len = $utgEnd - $utgBgn;
    die if ($len <= 0);

    my $per = $len / $length{$utg};

    if ($per < 0.95) {
        if ($length{$utg} < $TOOSHORT) {
            $spurShort++;
        } else {
            print STDERR "SPUR $utg len $length{$utg} mapped $len ($per%) coords $nucmer{$utg}\n";
            $spur++;
        }
    } else {
        $confirmed++;
        $confirmedBP += $length{$utg};
    }

    delete $nucmer{$utg};
}

print STDERR "FOUND ", scalar(keys %nucmer), " MATCHES, of which $spur are SPURs (and $spurShort are SHORT SPURS)\n";
print STDERR "FOUND $confirmed confirmed unitigs of total length $confirmedBP bp.\n";

########################################
#  Things with multiple matches are candidates for chimera.

my $classFalseMateSplit = 0;       my $lenFalseMateSplit = 0;
my $classShort = 0;                my $lenShort = 0;
my $classDisconnected = 0;         my $lenDisconnected = 0;
my $classRepeat = 0;               my $lenRepeat = 0;
my $classRepeatEndWithRepeat = 0;  my $lenRepeatEndWithRepeat = 0;
my $classRepeat1Unique = 0;        my $lenRepeat1Unique = 0;
my $classRepeat2Unique = 0;        my $lenRepeat2Unique = 0;
my $classMultipleBlock = 0;        my $lenMultipleBlock = 0;
my $classUnknown = 0;              my $lenUnknown = 0;


foreach my $utg (keys %nucmer) {
    my @m = split '\n', $nucmer{$utg};

    my $numBlocks = 0;
    my @bgn;
    my @end;

    my $min = 999999999;
    my $max = 0;

    #print STDOUT "$utg\n";

    #  Examine the matches, decide if the multiple matches are due to
    #    unitig is a repeat
    #    short repeats interior to the unitig
    #    long repeats at the end
    #    chimer

    my $numMatches = 0;
    my $isDisconnected = 0;
    my $verboseBuild = 1;

    print STDOUT "--------------------------------------------------------------------------------\n" if ($verboseBuild);
    print STDOUT "UNITIG $utg\n" if ($verboseBuild);

    foreach my $m (@m) {
        my ($utgBgn, $utgEnd, $genBgn, $genEnd) = split '\s+', $m;

        $numMatches++;
        print STDOUT "   $m" if ($verboseBuild);

        #  Are we an exact match to the largest thing so far?

        if (($utgBgn == $min) && ($utgEnd == $max)) {
            print STDOUT " -- is exact\n" if ($verboseBuild);
            next;
        }

        #  Search for things comtained or overlapping previous matches.  This assumes that matches
        #  are sorted by begin coord.  We need to handle:
        #
        #  ----------  and  -----------     not    --------  or     -----
        #     ----               ---------       ------           ----------
        #
        #  and, the ones we're more-or-less really looking for:
        #
        #  --------       OR  --------
        #       --------               -------

        #  This isolates out the first case above.  We don't want to save it as a region,
        #  because it is completely contained in a previous region.
        if ($utgEnd <= $max) {
            print STDOUT " -- is contained in a region\n" if ($verboseBuild);
            next;
        }

        #  If we aren't even intersecting, we've found the third case.
        my $dlabel = "";
        if (($max > 0) && ($max < $utgBgn)) {
            my $dist =  $utgBgn - $max;
            $dlabel = " DISCONNECT $dist";
            $isDisconnected++
        }

        #  Finally, update the min/max extents.  This must be last.
        if ($max == 0) {
            $min = $utgBgn;
            $max = $utgEnd;
        }
        if (($utgBgn - 5 <= $min) && ($max < $utgEnd)) {  #  max STRICTLY LESS than utgEnd!
            pop  @end;
            push @end, $utgEnd;
            $max = $utgEnd;
            print STDOUT " -- extends a region\n" if ($verboseBuild);
            next;
        }
        if ($utgBgn < $min)   { $min = $utgBgn; }
        if ($max < $utgEnd)   { $max = $utgEnd; }

        push @bgn, $utgBgn;
        push @end, $utgEnd;

        $numBlocks++;

        print STDOUT " -- makes a new region$dlabel\n" if ($verboseBuild);
    }

    #  One more pass through to count repeats.  Above we just found the min/max extent of the
    #  alignments on the unitig.  Here we can count if alignments are interior to the extent, at the
    #  end, or completely spanning it.
    #
    #  This fails by counting interior disconnected chimer (e.g., -- -- --) as 'isInteriorRepeat'

    my $isMaxExact = 0;
    my $isInteriorRepeat = 0;
    my $isTerminalRepeat = 0;


    foreach my $m (@m) {
        my ($utgBgn, $utgEnd, $genBgn, $genEnd) = split '\s+', $m;

        #  New match covers the whole region.
        if (($utgBgn - 5 <= $min) &&
            ($max <= $utgEnd + 5)) {
            $isMaxExact++;
            next;
        }

        #  New match completely contained, near the extent.
        if (($utgBgn - 5 <= $min) ||
            ($max <= $utgEnd + 5)) {
            $isTerminalRepeat++;
            next;
        }

        #  Otherwise, the match must be interior.
        $isInteriorRepeat++;
    }



    print STDERR "UNDEF LENGTH $utg\n" if (!defined($length{$utg}));
    print STDERR "ZERO LENGTH $utg\n" if ($length{$utg} == 0);
    print STDERR "SHORT LENGTH $utg\n" if ($length{$utg} < 64);

    if (($numMatches == 2) &&
        ($isDisconnected > 0) &&
        ($isTerminalRepeat == 2)) {
        $classFalseMateSplit++;
        $lenFalseMateSplit += $length{$utg};
        print "$utg (len $length{$utg}) is a FALSE MATE SPLIT.\n" if ($verboseBuild);
        next;
    }

    if ($isDisconnected > 0) {
        $classDisconnected++;
        $lenDisconnected += $length{$utg};
        print "$utg (len $length{$utg}) is DISCONNECTED.\n" if ($verboseBuild);
        next;
    }

    if ($isMaxExact >= 2) {
        $classRepeat++;
        $lenRepeat += $length{$utg};
        print "$utg (len $length{$utg}) is a REPEAT.\n" if ($verboseBuild);
        next;
    }

    if ($max < $TOOSHORT) {
        $classShort++;
        $lenShort += $length{$utg};
        next;
    }

    if (($isMaxExact == 1) &&
        ($isTerminalRepeat > 0) &&
        ($isInteriorRepeat > 0)) {
        $classRepeatEndWithRepeat++;
        $lenRepeatEndWithRepeat += $length{$utg};
        if (0) {
            print "----------\n";
            print "$utg extent $min $max numBlocks $numBlocks numMatches $numMatches maxExact $isMaxExact interiorRepeat $isInteriorRepeat terminalRepeat $isTerminalRepeat disconnect $isDisconnected\n";
            print "$nucmer{$utg}\n";
            for (my $iii=0; $iii<scalar(@bgn); $iii++) {
                print "$iii -- $bgn[$iii] $end[$iii]\n";
            }
            print "$utg (len $length{$utg}) is REPEAT on END OF UNIQUE with a repeat inside.\n";
        }
        #die if ($numBlocks > 1);
        next;
    }

    if (($isMaxExact == 1) &&
        ($isTerminalRepeat == $numMatches - 1)) {
        $classRepeat1Unique++;
        $lenRepeat1Unique += $length{$utg};
        if (0) {
            print "----------\n";
            print "$utg extent $min $max numBlocks $numBlocks numMatches $numMatches maxExact $isMaxExact interiorRepeat $isInteriorRepeat terminalRepeat $isTerminalRepeat disconnect $isDisconnected\n";
            print "$nucmer{$utg}\n";
            for (my $iii=0; $iii<scalar(@bgn); $iii++) {
                print "$iii -- $bgn[$iii] $end[$iii]\n";
            }
            print "$utg (len $length{$utg}) is REPEAT WITH A BIT OF UNIQUE ON ONE END.\n";
        }
        next;
    }

    if (($isMaxExact == 1) &&
        ($isTerminalRepeat == 0) &&
        ($isInteriorRepeat == $numMatches - 1)) {
        $classRepeat2Unique++;
        $lenRepeat2Unique += $length{$utg};
        #print "$utg (len $length{$utg}) is REPEAT WITH A BIT OF UNIQUE ON BOTH ENDS.\n";
        next;
    }

    if ($numBlocks > 1) {
        $classMultipleBlock++;
        $lenMultipleBlock += $length{$utg};
        if (0) {
            print "----------\n";
            print "$utg extent $min $max numBlocks $numBlocks numMatches $numMatches maxExact $isMaxExact interiorRepeat $isInteriorRepeat terminalRepeat $isTerminalRepeat disconnect $isDisconnected\n";
            print "$nucmer{$utg}\n";
            for (my $iii=0; $iii<scalar(@bgn); $iii++) {
                print "$iii -- $bgn[$iii] $end[$iii]\n";
            }
            print "$utg (len $length{$utg}) is MULTIPLE BLCOKS.\n";
        }
        next;
    }

    $classUnknown++;
    $lenUnknown += $length{$utg};

    print "$utg (len $length{$utg}) is UNCLASSIFIED.\n";
    print "$utg extent $min $max numBlocks $numBlocks numMatches $numMatches maxExact $isMaxExact interiorRepeat $isInteriorRepeat terminalRepeat $isTerminalRepeat disconnect $isDisconnected\n";
    print "$nucmer{$utg}\n";

    #  Try to make some sense out of it.  ....
}

print STDERR "classFalseMateSplit        $classFalseMateSplit ($lenFalseMateSplit bp)\n";
print STDERR "classDisconnected          $classDisconnected ($lenDisconnected bp)\n";
print STDERR "classRepeat                $classRepeat ($lenRepeat bp)\n";
print STDERR "classShort                 $classShort ($lenShort bp)\n";
print STDERR "classRepeatEndWithRepeat   $classRepeatEndWithRepeat ($lenRepeatEndWithRepeat bp)\n";
print STDERR "classRepeat1Unique         $classRepeat1Unique ($lenRepeat1Unique bp)\n";
print STDERR "classRepeat2Unique         $classRepeat2Unique ($lenRepeat2Unique bp)\n";
print STDERR "classMultipleBlock         $classMultipleBlock ($lenMultipleBlock bp)\n";
print STDERR "classUnknown               $classUnknown ($lenUnknown bp)\n";
