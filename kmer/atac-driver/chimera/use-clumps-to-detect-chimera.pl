#!/usr/bin/perl

#  Takes a path to a properly formatted atac file, uses that file to
#  detect potential chimeric scaffolds.

use strict;

my $atacFile   = undef;
my $reference  = "A";
my $noiseLevel = 1;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg eq "-A") {
        $reference = "A";
    } elsif ($arg eq "-B") {
        $reference = "B";
    } elsif ($arg eq "-n") {
        $noiseLevel = shift @ARGV;
    } elsif (-e $arg) {
        $atacFile = $arg;
    } else {
        print STDERR "Unknown option (or input file) '$arg'\n";
    }
}

if (! -e "$atacFile") {
    print STDERR "usage: $0 [-A | -B] file.atac\n";
    print STDERR "  -A    use the first assembly as the reference (default)\n";
    print STDERR "  -B    use the second assembly as the reference\n";
    exit(1);
}

open(ATAC, "< $atacFile") or die;
my @ATAC = <ATAC>;
chomp @ATAC;
close(ATAC);

my %ATACtoUID;

foreach my $line (@ATAC) {
    if ($line =~ m/assemblyFile(\d)=(.*)$/) {
        chomp $line;

        my $sequenceFile;

        if (($1 == 1) && ($reference eq "B")) {
            $sequenceFile = $2;
        }
        if (($1 == 2) && ($reference eq "A")) {
            $sequenceFile = $2;
        }

        #  If not defined, we don't need to read in these ID's.

        if (defined($sequenceFile)) {
            $sequenceFile =~ s/.fasta/.info/;

            die "Failed to find info on '$sequenceFile'\n" if (! -e $sequenceFile);

            print STDERR "Reading ATAC to UID map for '$sequenceFile'\n";

            open(F, "< $sequenceFile");
            while (<F>) {
                if (m/^G/) {
                    my @vals = split '\s+', $_;
                    $ATACtoUID{$vals[2]} = $vals[13];
                }
            }
            close(F);
        }
    }
}


#  0 1 2           3  4       5     6  7 8            9   10 11 12 13
#  M u H4467431a11 r1 B35LC:0 56097 66 1 HUREF4:36734 812 66 -1 # 10867
#  M u H4467431a10 r1 B35LC:0 56163 29 1 HUREF4:36734 782 29 -1 # 10867

#  Note that our match below does not match the non-clump marker "-1"


#  Find the scaffolds with errors
#
#  Save the clump id for the first instance of every scaffold.  If
#  we've seen the scaffold before, and the clump id is now different,
#  remember this scaffold.

my %scaffold;
my %errors;

foreach (@ATAC) {
    if (m/^M\su\s.*\s#\s(\d+)$/) {
        my @v = split '\s+', $_;

        if (!defined($scaffold{$v[8]})) {
            $scaffold{$v[8]} = $v[13];
        } elsif ($scaffold{$v[8]} ne $v[13]) {
            $errors{$v[8]}++;
        }
    }
}

#  Print them
#
#  Go through the map again, remembering the number of times we see a
#  scaffold/clump pair.  It's also useful to remember the sum of the
#  lengths for this pair, and the chromosome it maps to.

my %counts;
my %length;
my %chrid;

foreach (@ATAC) {
    if (m/^M\su\s.*\s#\s(\d+)$/) {
        my @v = split '\s+', $_;

        if (defined($errors{$v[8]})) {
            my $string = "$v[8]\t$ATACtoUID{$v[8]}\t$v[13]";
            $counts{$string}++;
            $length{$string} += $v[10];
            $chrid{$string} = $v[4];
        }
    }
}

#  We could provide a raw dump of this data, but we'd like
#  to first denoise it.  A very simple denoser works - just
#  don't report anything with one match.

open(F, "| sort -k3,3");
foreach my $s (keys %counts) {
    if ($counts{$s} > $noiseLevel) {
        print F "$counts{$s}\t$length{$s}\t$s\t$chrid{$s}\n";
    }
}
close(F);
