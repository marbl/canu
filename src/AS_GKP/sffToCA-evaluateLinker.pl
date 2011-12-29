#!/usr/bin/perl
#
#  Run sffToCA on some sff file, then search the resulting reads for linker.
#

###
###  Cobbled together from bits and pieces leftover from a one-off analysis.
###                METHOD WORKS, SCRIPT PROBABLY WILL NOT.
###

use strict;

my $sffFile       = shift @ARGV;
my $linkerType    = shift @ARGV;
my $outputPrefix  = shift @ARGV;

#  Insert size and library name don't matter.

if (! -e "$outputPrefix.frg") {
    my $cmd;
    $cmd  = "sffToCA";
    $cmd .= " -insertsize 3000 300 -libraryname T -clear 454 -trim chop -linker $linkerType";
    $cmd .= " -output $outputPrefix";
    $cmd .= " $sffFile";
    system($cmd);
}

#  This is more or less AS_REF/frg-to-fasta.pl

if (! -e "$outputPrefix.fasta") {
    my $acc;

    open(F, "< $outputPrefix.frg")   or die "Failed to open '$outputPrefix.frg' for converting to fasta.\n";
    open(O, "< $outputPrefix.fasta") or die "Failed to open '$outputPrefix.fasta' for writing.\n";

    while (<F>) {
        chomp;

        if (m/acc:(.*)$/) {
            $acc = $1;
        }
        if (m/seq:/) {
            print O ">$acc\n";
            $_ = <F>;
            chomp;

            while (!m/^.$/) {
                print O "$_";
                $_ = <F>;
                chomp;
            }

            print O "\n";
        }
    }

    close(O);
    close(F);
}

####
####  Align those fragments to reference.
####

if (! -e "$outputPrefix.sim4db") {
}

if (! -e "$outputPrefix.sim4db.extent") {
}

if (! -e "") {
}

####
####
####  This was called 'extract-unaligned.pl'.  It seems to mask any sequence that
####  aligns with N's, then extracts the unaligned bases.
####
####

#!/usr/bin/perl

use strict;

my %seq;

open(F, "< FRAG.fasta") or die;
while (!eof(F)) {
    my $hdr = <F>;  chomp $hdr;  $hdr =~ s/>//;
    my $seq = <F>;  chomp $seq;

    $seq{$hdr} = $seq;
    $len{$hdr} = length($seq);
}
close(F);


open(M, "convertToExtent < FRAG-ms16.sim4db |") or die;
$_ = <M>;  #  Header line.
while (<M>) {
    my @v = split '\s+', $_;

    my $seq = $seq{$v[0]};
    die if (!defined($seq));

    my $bgn = ($v[3] < $v[4]) ? $v[3] : $v[4];
    my $end = ($v[3] < $v[4]) ? $v[4] : $v[3];

    #print  "$bgn $end\n";
    #print  "$seq\n";
    substr($seq, $bgn, $end - $bgn) = "N" x ($end - $bgn);
    #print  "$seq\n";

    $seq{$v[0]} = $seq;
}
close(M);

foreach my $k (keys %seq) {
    my $s = $seq{$k};
    my @s = split 'N+', $s;
    my $i = "00";

    foreach my $S (@s) {
        if (length($S) > 6) {
            print ">$k-$i\n$S\n";
            $i++;
        }
    }
}


####
####
####  OLDER method, using nucmer.  Too slow.
####
####

#
#  Align fragments to the reference with nucmer.  WARNING!  The 'postnuc' step takes about TWO HOuRS
#  on a half-plate of FLX mates (even EXCLUDING the mated reads, just the fragment reads).
#

if (! -e "$prefix.nucmer") {
    my $cmd;

    #  extend delta optimize simplify
    #  --delta    - defalut - generate .delta file (DO_DELTA = true) (FAILS if --nodelta used)
    #  --extend   - default - cluster extension (DO_EXTEND = true)
    #  --optimize - default - align score optimization - if hit end of sequence, reverse dp to optimize score from there (DO_SEQEND = false)
    #  --simplify - default - remove shadowed clusters (alignments to itself?) (DO_SHADOWS = false)
    #
    $cmd  = "nucmer";
    $cmd .= " --maxmatch";
    $cmd .= " --coords";
    $cmd .= " -p $prefix";
    $cmd .= " $prefix.fasta";
    $cmd .= " $refName";

    system($cmd);
}

#
#  Parse alignments
#
#
# start1     end1        BAR       start2    end2      BAR       len1      len2      BAR       identity  BAR       refID     qryID     

my $perfect = 0;
my $partial = 0;

my %seen;
my %cnts;

my %fulll;
my %lpart;
my %mpart;
my %rpart;

open(F, "< $prefix.coords") or die "Failed to open '$prefix.coords' for reading.\n";
$_ = <F>;
$_ = <F>;
$_ = <F>;
$_ = <F>;
$_ = <F>;
while (<F>) {
    tr/^\s+//;
    tr/\s+$//;

    my @v = split '\s+', $_;

    my $bgn = ($v[4] < $v[5]) ? $v[4] : $v[5];
    my $end = ($v[4] < $v[5]) ? $v[5] : $v[4];

    $bgn--;
    $end = $frgLen{$v[13]} - $end;

    if (($bgn == 0) &&
        ($end == 0)) {
        $perfect++;
    } else {
        $partial++;
    }

    $seen{$v[13]} .= "$_";
    $cnts{$v[13]}++;

    if      (($bgn < 40) && ($end < 40)) {
        #  Full length
        $fulll{$v[13]}++;

    } elsif ($bgn < 40) {
        #  Left side
        $lpart{$v[13]}++;

    } elsif ($end < 40) {
        #  Right side
        $rpart{$v[13]}++;

    } else {
        #  No side
        $mpart{$v[13]}++;
    }
}
close(F);

foreach my $k (keys %cnts) {
    if (($cnts{$k} > 1) &&
        ($fulll{$k} == 0) &&
        ($mpart{$k} == 0)) {
        print "----------------------------------------\n";
        print "$k -- len $frgLen{$k} -- $lpart{$k} -- $rpart{$k}\n";
        print "$seen{$k}";
    }
}

print "PERFECT $perfect\n";
print "PARTIAL $partial\n";
