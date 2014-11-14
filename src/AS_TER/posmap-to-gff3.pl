#!/usr/bin/perl

use strict;

my $dir = "asm/9-terminator";
my $asm = "asm";

my %scfOfInterest;
my %ctgOfInterest;
my %utgOfInterest;

$scfOfInterest{"scf7180002220303"}++;
$scfOfInterest{"scf7180002224952"}++;
$scfOfInterest{"scf7180002258291"}++;
$scfOfInterest{"scf7180002259533"}++;
$scfOfInterest{"scf7180002280184"}++;
$scfOfInterest{"scf7180002290117"}++;
$scfOfInterest{"scf7180002310267"}++;
$scfOfInterest{"scf7180002311635"}++;

print "#gff-version 3\n";

open(F, "< $dir/$asm.posmap.scflen") or die;
while (<F>) {
    my @v = split '\s+', $_;
    if ((scalar(keys %scfOfInterest) > 0) && (exists($scfOfInterest{"scf$v[0]"}))) {
        printf "scf$v[0]\t$asm\tca_scaffold\t1\t$v[1]\t.\t.\t.\tID=scf$v[0];Name=scf$v[0]\n";
    }
}
close(F);


open(F, "< $dir/$asm.posmap.ctgscf") or die;
while (<F>) {
    my @v = split '\s+', $_;
    if ((scalar(keys %scfOfInterest) > 0) && (exists($scfOfInterest{"scf$v[1]"}))) {
        $v[2]++;
        $v[4] = "+" if ($v[4] eq "f");
        $v[4] = "-" if ($v[4] eq "r");

        print "scf$v[1]\t$asm\tca_contig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tID=ctg$v[0];Name=ctg$v[0];Parent=scf$v[1]\n";

        $ctgOfInterest{"ctg$v[0]"}++;
    }
}
close(F);

#  LOAD UNITIG MARKINGS.  This is a temporary hack until posmap is updated.
my %utgStatus;
open(F, "< $dir/$asm.utgstatus") or die;
while (<F>) {
    my @v = split '\s+', $_;
    $utgStatus{$v[1]} = $v[8];
}
close(F);


open(F, "< $dir/$asm.posmap.utgscf") or die;
while (<F>) {
    my @v = split '\s+', $_;
    if ((scalar(keys %scfOfInterest) > 0) && (exists($scfOfInterest{"scf$v[1]"}))) {
        $v[2]++;
        $v[4] = "+" if ($v[4] eq "f");
        $v[4] = "-" if ($v[4] eq "r");

        #print "scf$v[1]\t$asm\tca_unitig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tID=utg$v[0];Name=utg$v[0]\n";

        if      ($utgStatus{$v[0]} eq "N") {
            #  Unplaced, degenerate
            #print "scf$v[1]\t$asm\tca_degenerate_unitig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tParent=utg$v[0]\n";
            print "scf$v[1]\t$asm\tca_degenerate_unitig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tID=deg$v[0]\n";
        } elsif ($utgStatus{$v[0]} eq "U") {
            #  Placed, unique
            #print "scf$v[1]\t$asm\tca_unique_unitig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tParent=utg$v[0]\n";
            print "scf$v[1]\t$asm\tca_unique_unitig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tID=deg$v[0]\n";
        } elsif ($utgStatus{$v[0]} eq "S") {
            #  Placed, surrogate
            #print "scf$v[1]\t$asm\tca_surrogate_unitig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tParent=utg$v[0]\n";
            print "scf$v[1]\t$asm\tca_surrogate_unitig\t$v[2]\t$v[3]\t.\t$v[4]\t.\tID=deg$v[0]\n";
        } else {
            die "Unknown unitig status $_\n";
        }
    }
}
close(F);


open(F, "< $dir/$asm.posmap.frgscf") or die;
while (<F>) {
    my @v = split '\s+', $_;
    if ((scalar(keys %scfOfInterest) > 0) && (exists($scfOfInterest{"scf$v[1]"}))) {
        $v[2]++;
        $v[4] = "+" if ($v[4] eq "f");
        $v[4] = "-" if ($v[4] eq "r");

        print "scf$v[1]\t$asm\tca_fragment\t$v[2]\t$v[3]\t.\t$v[4]\t.\tID=frg$v[0];Name=frg$v[0]\n";
    }
}
close(F);


open(F, "< $dir/$asm.posmap.varscf") or die;
while (<F>) {
    my @v = split '\s+', $_;
    if ((scalar(keys %scfOfInterest) > 0) && (exists($scfOfInterest{"scf$v[1]"}))) {
        print "scf$v[1]\t$asm\tca_variant\t$v[2]\t$v[3]\t.\t.\t.\n";
    }
}
close(F);
