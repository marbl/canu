#!/usr/bin/perl

#  examines an atac mapping, counts the number of times a scaffold is
#  mapped to wildly different places -- wildly being more than a few
#  bp away on the same chromosome (mind the gap, please!) or (gasp!) a
#  different chromosome.

use strict;

#  scafA is really scaffolds-from-map1, and scafB is scaffolds-from-map2.

my %scafA_to_chr;
my %scafB_to_chr;

my %scafA_to_chr_mismatch;
my %scafB_to_chr_mismatch;

my $scafAlen;
my $scafBlen;

open(F, "< overlap.Aannotation");
while (<F>) {
    chomp;

    if (m/^N\s+(\d+):(\d+)-(\d+)\[\s*\d+\].*\s(\d+):\s*(\d+)-\s*(\d+)\).*\s(\d+):\s*(\d+)-\s*(\d+)\)/) {
        my $id1 = $1;
        my $b1  = $2;
        my $e1  = $3;

        my $id2a = $4;
        my $b2a  = $5;
        my $e2a  = $6;

        my $id2b = $7;
        my $b2b  = $8;
        my $e2b  = $9;


        if (($id2a > 0) && ($e2a > 0)) {
            if (defined($scafA_to_chr{$id2a})) {
                if ($scafA_to_chr{$id2a} != $id1) {
                    if (! defined($scafA_to_chr_mismatch{$id2a})) {
                        $scafA_to_chr_mismatch{$id2a} = $scafA_to_chr{$id2a};
                    }
                    $scafA_to_chr_mismatch{$id2a} .= "\1$id1\0$_";
                }
            } else {
                $scafA_to_chr{$id2a} = "$id1\0$_";
            }
        }

        if (($id2b > 0) && ($e2b > 0)) {
            if (defined($scafB_to_chr{$id2b})) {
                if ($scafB_to_chr{$id2b} != $id1) {
                    if (! defined($scafB_to_chr_mismatch{$id2b})) {
                        $scafB_to_chr_mismatch{$id2b} = $scafB_to_chr{$id2b};
                    }
                    $scafB_to_chr_mismatch{$id2b} .= "\1$id1\0$_";
                }
            } else {
                $scafB_to_chr{$id2b} = "$id1\0$_";
            }
        }
    }
}
close(F);

#  Count the number of things in *_mismatch that are the same
#
my %merge;
my $both;
foreach my $f (keys %scafA_to_chr_mismatch) {
    $merge{$f}++;
}
foreach my $f (keys %scafB_to_chr_mismatch) {
    $merge{$f}++;
}
foreach my $f (keys %merge) {
    if (defined($scafA_to_chr_mismatch{$f}) && defined($scafB_to_chr_mismatch{$f})) {
        $both++;
    }
}

print "scafA count: ", scalar(keys %scafA_to_chr_mismatch), "\n";
print "scafB count: ", scalar(keys %scafB_to_chr_mismatch), "\n";
print "both  count: ", $both, "\n";


#  Run through the input again, pulling out matches that map a single
#  scaffold to two different chromosomes, then what?

if (0) {

#  flawed method for doing this

open(F, "| sort -k6n -k7n > scaffold-consistency.map1dups");
foreach my $f (values %scafA_to_chr_mismatch) {
    my @things = split "\1", $f;
    foreach my $t (@things) {
        my ($a, $b) = split '\0', $t;
        print F "$b\n";
    }
}
close(F);

open(F, "| sort -k11n -k12n > scaffold-consistency.map2dups");
foreach my $f (values %scafB_to_chr_mismatch) {
    my @things = split "\1", $f;
    foreach my $t (@things) {
        my ($a, $b) = split "\0", $t;
        print F "$b\n";
    }
}
close(F);

}

#  we saved the iid of the scaffold that maps to different chromosomes as the key,
#  so just parse the matches again, pulling out all those scaffolds.

open(A, "| sort -k6n -k7n > scaffold-consistency.map1dups");
open(B, "| sort -k11n -k12n > scaffold-consistency.map2dups");
my $matchesA = 0;
my $matchesB = 0;
open(F, "< overlap.Aannotation");
while (<F>) {
    if (m/^.\s+(\d+):(\d+)-(\d+)\[\s*\d+\].*\s(\d+):\s*(\d+)-\s*(\d+)\).*\s(\d+):\s*(\d+)-\s*(\d+)\)/) {
        my $id1 = $1;
        my $b1  = $2;
        my $e1  = $3;

        my $id2a = $4;
        my $b2a  = $5;
        my $e2a  = $6;

        my $id2b = $7;
        my $b2b  = $8;
        my $e2b  = $9;

        if (defined($scafA_to_chr_mismatch{$id2a}) || defined($scafA_to_chr_mismatch{$id2b})) {
            print A $_;
            $matchesA++;
        }

        if (defined($scafB_to_chr_mismatch{$id2a}) || defined($scafB_to_chr_mismatch{$id2b})) {
            print B $_;
            $matchesB++;
        }
    }
}

print STDERR "matches for map1: $matchesA\n";
print STDERR "matches for map2: $matchesB\n";
