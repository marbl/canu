#!/usr/bin/perl

use strict;

#  Examines an atac mapping, counts the number of times a scaffold is
#  mapped to wildly different places -- wildly being more than a few
#  bp away on the same chromosome (mind the gap, please!) or (gasp!) a
#  different chromosome.
#
#  Assumes that the Aannotation primary axis is chromosomes.
#
#  Change the first . in the m// below to restrict to specific types
#  of regions, e.g., N's.  Useful choices here are:
#    .    - all regions
#    U    - unmapped (will do nothing)
#    1    - has only one mapping
#    Y    - just those they agree on
#    N    - disagree, but on the same destination
#    !    - disagree, and on different destinations
#    ?    - inconsistent mapping (OK, this one isn't useful)
#
#    1Y   - the regions that have no disagreement
#    1YN! - all consistent regions

if (scalar(@ARGV != 2)) {
    print STDERR "usage: $0 some.Aannotation outprefix\n";
    exit(1);
}

my $filename   = shift @ARGV;
my $outprefix  = shift @ARGV;

#  scafA is really scaffolds-from-map1, and scafB is scaffolds-from-map2.

my (%scafA_to_chr, %scafA_to_chr_mismatch, $scafAlen);
my (%scafB_to_chr, %scafB_to_chr_mismatch, $scafBlen);

open(F, "< $filename");
while (<F>) {
    chomp;

    if (m/^[1YN!]\s+(\d+):(\d+)-(\d+)\[\s*\d+\].*\s(\d+):\s*(\d+)-\s*(\d+)\).*\s(\d+):\s*(\d+)-\s*(\d+)\)/) {
        my ($id1,  $b1,  $e1)  = ($1, $2, $3);
        my ($id2a, $b2a, $e2a) = ($4, $5, $6);
        my ($id2b, $b2b, $e2b) = ($7, $8, $9);

        #  If we have a mapping from method A or method B, save the
        #  chromosome that the scaffold mapped to.  If we've already
        #  mapped this scaffold to some other chromosome, call it a
        #  mismatch.

        if (($id2a > 0) && ($e2a > 0)) {
            if (defined($scafA_to_chr{$id2a})) {
                if ($scafA_to_chr{$id2a} != $id1) {
                    $scafA_to_chr_mismatch{$id2a}  = $scafA_to_chr{$id2a} if (! defined($scafA_to_chr_mismatch{$id2a}));
                    $scafA_to_chr_mismatch{$id2a} .= "\1$id1\0$_";
                }
            } else {
                $scafA_to_chr{$id2a} = "$id1\0$_";
            }
        }

        if (($id2b > 0) && ($e2b > 0)) {
            if (defined($scafB_to_chr{$id2b})) {
                if ($scafB_to_chr{$id2b} != $id1) {
                    $scafB_to_chr_mismatch{$id2b}  = $scafB_to_chr{$id2b} if (! defined($scafB_to_chr_mismatch{$id2b}));
                    $scafB_to_chr_mismatch{$id2b} .= "\1$id1\0$_";
                }
            } else {
                $scafB_to_chr{$id2b} = "$id1\0$_";
            }
        }
    }
}
close(F);

#  Count the number of things in *_mismatch that are the same.
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

print "num scafA: ", scalar(keys %scafA_to_chr_mismatch), "\n";
print "num scafB: ", scalar(keys %scafB_to_chr_mismatch), "\n";
print "num both:  ", $both, "\n";

#  Run through the input again, pulling out matches that map a
#  single scaffold to two different chromosomes, then what?  we
#  saved the iid of the scaffold that maps to different
#  chromosomes as the key, so just parse the matches again,
#  pulling out all those scaffolds.

open(A, "| sort -k6n  -k7n > $outprefix.scaffold-consistency.map1dups");
open(B, "| sort -k11n -k12n > $outprefix.scaffold-consistency.map2dups");
my $matchesA = 0;
my $matchesB = 0;
open(F, "< $filename");
while (<F>) {
    if (m/^.\s+(\d+):(\d+)-(\d+)\[\s*\d+\].*\s(\d+):\s*(\d+)-\s*(\d+)\).*\s(\d+):\s*(\d+)-\s*(\d+)\)/) {
        my ($id1,  $b1,  $e1)  = ($1, $2, $3);
        my ($id2a, $b2a, $e2a) = ($4, $5, $6);
        my ($id2b, $b2b, $e2b) = ($7, $8, $9);

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
close(A);
close(B);

print STDERR "matches for map1: $matchesA\n";
print STDERR "matches for map2: $matchesB\n";
