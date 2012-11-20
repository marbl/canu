#!/usr/bin/perl

use strict;

#  Simple hack to merge multiple QC reports into one.  Uses the first
#  QC on the command line to figure out what fields exist, then writes
#  a multiple column QC using those fields and data from all QC files
#  on the command line.

#  BUG: if one file is empty (say, the second file) then NO values
#  are added to the columns, and we have just shifted all the values.
#
#  In this example, MU1/MU1.qc is actually the empty one, but it looks
#  like OU1/OU1.qc is it.
#
#  [files]
#                          MB1/MB1.qc      MU1/MU1.qc      OB1/OB1.qc      OU1/OU1.qc
#
#  [Scaffolds]
#  TotalScaffolds          1954            2046            3459
#  TotalContigsInScaffo    1954            2046            3459
#  MeanContigsPerScaffo    1.00            1.00            1.00
#  MinContigsPerScaffol    1               1               1

my @labels;
my %values;

my $s;
my $l;
my $v;

my $firstFile = 1;
my $numFiles  = 1;  #  First file is the label column
my $isWiki    = 0;
my $isCSV     = 0;

push @labels, "Files";

if ($ARGV[0] eq "-wiki") {
    shift @ARGV;

    $isWiki = 1;
    print "{| class=\"wikitable\" border=\"1\"\n";
}

if ($ARGV[0] eq "-csv") {
    shift @ARGV;

    $isCSV = 1;
}

while (scalar(@ARGV) > 0) {
    if (open(F, "< $ARGV[0]")) {

        if ($isWiki) {
            $values{"Files"} .= "|| $ARGV[0]";
            $values{"Files"} .= "BRIWASHERE";
        } elsif ($isCSV) {
            $values{"Files"} .= ",$ARGV[0]";
            $values{"Files"} .= "BRIWASHERE";
        } else {
            $values{"Files"} .= "\t$ARGV[0]";
            $values{"Files"} .= "BRIWASHERE";
        }

        $numFiles++;

        while (<F>) {
            $_ =~ s/^\s+//;
            $_ =~ s/\s+$//;

            if ($_ =~ m/^\s+$/) {
                next;
            }

            if ($_ =~ m/^\[(.*)\]$/) {
                $s = $1;
                next;
            }

            if ($_ =~ m/^(.*)=(.*)$/) {
                $l = "$s\0$1";
                $v = $2;
            } else {
                next;
            }

            if ($isCSV) {
                $l =~ s/,/ /g;
                $v =~ s/,/ /g;
            }

            if ($firstFile) {
                push @labels, $l;
            }

            if ($isWiki) {
                $values{$l} .= "|| $v";
                $values{$l} .= "BRIWASHERE";
            } elsif ($isCSV) {
                $values{$l} .= ",$v";
                $values{$l} .= "BRIWASHERE";
            } else {
                $values{$l} .= substr("\t$v                ", 0, 16);
                $values{$l} .= "BRIWASHERE";
            }
        }
        close(F);

        my @k = keys %values;
        foreach my $l (@k) {
            if ($values{$l} =~ m/^(.*)BRIWASHERE$/) {
                $values{$l}  = $1;
            } else {
                if ($isWiki) {
                    $values{$l} .= "|| N/A";
                } elsif ($isCSV) {
                    $values{$l} .= ",0";
                } else {
                    $values{$l} .= substr("\tN/A                ", 0, 16);
                }
            }
        }

        $firstFile = 0;
    }

    shift @ARGV;
}


my $lasts;
foreach my $xx (@labels) {
    ($s, $l) = split '\0', $xx;

    if ($s ne $lasts) {
        if ($isWiki) {
            print "| colspan=\"$numFiles\" rowheight=\"6\" bgcolor=\"#f2f2f2\" |\n";
            print "|-\n";
            print "| colspan=\"$numFiles\" bgcolor=\"#d2d2d2\" | $s\n";
            print "|-\n";
        } else {
            print "\n[$s]\n";
        }
        $lasts = $s;
    }

    if ($isWiki) {
        print "| $l $values{$xx}\n";
        print "|-\n";
    } elsif ($isCSV) {
        print "$l$values{$xx}\n";
    } else {
        $l = substr("$l                    ", 0, 20);
        print "$l$values{$xx}\n";
    }
}

if ($isWiki) {
    print "|}\n";
}
