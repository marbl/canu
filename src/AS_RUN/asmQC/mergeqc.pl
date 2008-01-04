#!/usr/bin/perl

#  Simple hack to merge multiple QC reports into one.  Uses the first
#  QC on the command line to figure out what fields exist, then writes
#  a multiple column QC using those fields and data from all QC files
#  on the command line.

use strict;

my @labels;
my %values;

my $s;
my $l;
my $v;

my $firstFile = 1;

push @labels, "files";

while (scalar(@ARGV) > 0) {
    open(F, "< $ARGV[0]") or die "Failed to open '$ARGV[0]'\n";

    $values{"files"} .= "\t$ARGV[0]";

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

        if ($firstFile) {
            push @labels, $l;
        }

        $values{$l} .= substr("\t$v                ", 0, 16);
    }
    close(F);

    $firstFile = 0;

    shift @ARGV;
}


my $lasts;
foreach my $xx (@labels) {
    ($s, $l) = split '\0', $xx;

    if ($s ne $lasts) {
        print "\n[$s]\n";
        $lasts = $s;
    }

    $l = substr("$l                    ", 0, 20);

    print "$l$values{$xx}\n";
}
