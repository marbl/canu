#!/usr/local/bin/perl

=head1 NAME

glen_est_truncadjusted.pl - estimate genome length based on overlap statistics.

=head1 VERSION

This document refers to version 1.00 of glen_est_truncadjusted.pl, released 05.01.2006.

=head1 SYNOPSIS

glen_est_truncadjusted.pl -base <asm prefix> -K <int> [-first <int>] [-last <int>] [-into <int>] [-window <int>] [-minovl <int>]

Flags:

    -base <asm base name> Path plus prefix of Celera Assembler {frg,gkp,ovl}Stores
    -K <int>              Exclude fragments with more than this many overlaps
    -first <int>          First fragment to analyze; default=1
    -last <int>           Last fragment to analyze; default=100000
    -into <int>           Length of skipped prefix of analyzed fragments; default=50
    -window <int>         Analysis window width; default=100
    -minovl <int>         Minimum overlap length; default=40
    -help|h               Print help message

Assumptions:
    No deleted fragments in fragStore
    K value at least twice the average depth of coverage in unique regions
    Values of first,last reasonable given size of dataset
    Value of into+window is less than length of most fragments
    Fragments in [first,last] are typical of entire dataset
    Single genome with little contamination

=head1 DESCRIPTION

=head2 Overview

Estimates genome length for a sequencing project from a low or
intermediate coverage dataset, discounting fragments whose coverage is
high enough to suggest they may be repetitive and avoiding near
duplicate fragments in case of artifact.

=head2 Credit

Initial version described by Art Delcher (Kirkness et al 2004,
Science).  Modifications and coding of current version by Aaron
Halpern.

=cut


use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);




my %options = ();
my $prog = $0;
$prog =~ s/.*\///;
my ($help,$man,$base,$ovlStore,$frgStore,$gkpStore,$K,$first,$last,$Into,$windowsize,$Minovl);
$first=1;
$last=100000;
$Into=50;
my $K=0;
my $windowsize=100;
my $Minovl=40;
my $caRoot = $AS_ROOT;
my $caBin = "";

BEGIN: {
  GetOptions(
    "help|h"      => \$help,
    "base=s"         => \$base,
    "K=i"     => \$K,
    "first=i" => \$first,
    "last=i" => \$last,
    "into=i" => \$Into,
    "window=i" => \$windowsize,
    "minovl=i" => \$Minovl,
    "CA=s" => \$caRoot,
    "CABIN=s" => \$caBin
  ) || pod2usage(2);
  pod2usage(1) if defined($help) || !defined($base) || $K<=0 || $first >= $last || $first < 1 || $Into < 0 || $Into > 500 || $windowsize<=0 || $windowsize>400;
}

# simple factorials
my @fact;
my $f;
$fact[0]=1;
$fact[1]=1;
for($f=2;$f<$K;$f++){
    $fact[$f]=$fact[$f-1]*$f;
}

#$frgStore = $base . ".frgStore";
#frg store no longer exists, part of gkp store
$gkpStore = $base . ".gkpStore";
$ovlStore = $base . ".ovlStore";

if ($caBin eq "") {
   if($caRoot eq ""){
       $caRoot="/bioinfo/work/software/released/asm/current";
   }

   if( $ENV{OSTYPE} eq "osf1" ){
   	$caBin="$caRoot/OSF1";
   }
   elsif ($ENV{OSTYPE} eq "linux"){
       if ($ENV{HOSTTYPE} eq "x86_64") {
         $caBin="$caRoot/Linux-amd64";
       }
       else {
          $caBin="$caRoot/Linux-i686";
       }
   }
   else {
      die "Unknown OSTYPE $OSTYPE\n";
   }
}

open(ALLFRG,"${caBin}/bin/gatekeeper -lastfragiid $gkpStore|");
$tmpstr=<ALLFRG>;
close(ALLFRG);
chomp $tmpstr;
my @foo = split / /,$tmpstr;
my $nfoo = @foo;
my $nfrags=$foo[$nfoo-1];

print STDERR "The last frag id is $nfrags\n";

if($first >$nfrags){
    print STDERR "First must be smaller than number of fragments!\n";
    pod2usage(1);
}

if($last >$nfrags){
    $last=$nfrags;
    print STDERR "WARNING: Truncating end of fragment range: Store not that big!\n";
}

my $avgfrglen=0;
my $nfr=0;
my $nzero=0;
open(LENS,"${caBin}/bin/gatekeeper -b $first -e $last -clear LATEST -dumpfragments -tabular $gkpStore | awk '{print \$1\"\\t\"\$13-\$12}' |");
print "Running ${caBin}/bin/gatekeeper -b $first -e $last -clear LATEST -dumpfragments -tabular $gkpStore | awk '{print \$1\"\\t\"\$13-\$12}'\n";
while(<LENS>){
    @w=split;

    if ($w[0] eq "UID") { next; }
    $avgfrglen+=$w[1];
    $nfr++;
    if($w[1]==0){
      $nzero++;
}
}
close(LENS);

if($nfr != $last-$first+1 || $nfr-$nzero == 0){
    print STDERR "Trouble with fragment interval!\n";
    exit(-1);
}
$nfr-=$nzero;
$avgfrglen/=($nfr-$nzero);

print "Running command ${caBin}/bin/overlapStore -b $first -e $last -d $ovlStore\n";
open(OLAPS,"${caBin}/bin/overlapStore -b $first -e $last -d $ovlStore -g $gkpStore|");
my $NR=0;
my $prev=-1;
my $readswithovls=0;
my $n=0;
my $nonrepeatreads=0;
my $nonrepeatovls=0;
my $N=0;

#my @ovlCount;
#my @ovlCountInWin;
#for($i=0;$i<$last;$i++){
#    $ovlCount[$i]=$n;
#    $ovlCountInWin[$i]=$N;
#}

while(<OLAPS>){
    my @w=split;
    if($w[3]<=$Into){next;}
    $NR++;
    if($w[0]!=$prev){
      $readswithovls++;
#      $ovlCount[$prev-$first]=$n;
#      $ovlCountInWin[$prev-$first]=$N;
      if($n<$K&&$NR>1){
	$nonrepeatreads++;
        $nonrepeatovls+=$N;
      }
      $n=0;
      $N=0;
    }
    $prev=$w[0];
    $n++;
    if($w[3]<=$Into+$windowsize){
	$N++;
    }
}
close(OLAPS);

if($prev>0){
#    $ovlCount[$prev-$first]=$n;
#    $ovlCountInWin[$prev-$first]=$N;
    if($n<$K){
	$nonrepeatreads++;
	$nonrepeatovls+=$N;
    }
}

if( $nonrepeatreads/$nfr < .7 ) {
    my $usedFrac=1-$nonrepeatreads/$nfr;
    print "SERIOUS WARNING! More than thirty percent ($usedFrac) of data marked repeat!\n";
    print "\tYou may have specified too small a value for K, which\n";
    print "\tcan lead to a BIG misestimation of genome length!\n";
} else {
    my $goodfrac= $nonrepeatreads/$nfr;
    print STDERR "Using $goodfrac of analyzed reads\n";
}

my $lambda=$nonrepeatovls/($nfr-$readswithovls+$nonrepeatreads);

my $adjusted_lambda=compute_truncated_mean($K,$lambda,$Into,$avgfrglen,$Minovl,$windowsize);

my $glen = $windowsize * $nfrags / $lambda;

my $fulllambda= $lambda * $avgfrglen/$windowsize;

print "final estimate $nfrags $fulllambda $glen\n";

exit(0);



### Need to adjust naive estimate of lambda to deal with the fact
### that if K is small enough, we have filtered out a significant
### amount of upper tail of unique coverage.  The following subroutine
### computes, by more or less Newton_s method, that poisson mean
### that, when truncated >= K, would give (truncated) mean == naive
### lambda estimate

sub compute_truncated_mean(){
    my $precision=.00001;

    # truncate at less than this many fragments overlapping the read
    my $k=shift;

    # goal value for fitted parameter
    my $goal=shift;

    # how far into genome?
    my $into=shift;

    # average fragment size
    my $frglen=shift;

    # minimum ovl thickness
    my $minovl=shift;

    #window width
    my $width=shift;

    # scale-factor for how much of read is evaluated for k and how much counts ovls
    my $frac=815;
    $frac-=$minovl;
    $frac-=$into;
    $frac/=$width;

    my $high=$goal*10;
    my $highest=estimate($high,$frac,$k);
    my $low=$precision;
    my $lowest=estimate($low,$frac,$k);
    my $mid=($high+$low)/2;
    my $midest=estimate($mid,$frac,$k);

    if($highest<$goal){
    print STDERR "Trouble getting initial bounds!\n";
    exit(-1);
    }
    my $err=abs($midest-$goal);
    while($err>$precision){
	if($midest>$goal){
	    $highest=$midest;
	    $high=$mid;
        } else {
	    $lowest=$midest;
	    $low=$mid;
	}
        $mid=($high+$low)/2;
        $midest=estimate($mid,$frac,$k);
        $err=abs($midest-$goal);
    }
    return ($mid);
}

### core: compute the mean overlap count after truncation given
### a poisson distrib with mean = $l
sub estimate(){

    my $l=shift;
    my $frac=shift;
    my $k=shift;

    my $lambda=$l*$frac;  # this is the true mean in the full fragment

    my $L=0;  # this is the estimated mean in the first 100 bp

    my $P = 0; # normalizing probability ... see below

    # compute mean by considering full fragment from 1 to k
    my $i;
    for($i=0;$i<$k;$i++){

	# for i fragments starting during the read

	# prob of i is ...
	my $prob=exp($i*log($lambda))*exp(-$lambda)/$fact[$i];

	# and expected number occuring within the first 100 bp is ...
	my $y=$i/$frac;

	# so the contribution to the expected value is ...
	$L+=$y*$prob;

	# but since we exclude the cases affected truncation, need
	# to compute the denominator
	$P+=$prob;
    };
    return($L/$P);
}

sub max(){
    $a=shift;
    $b=shift;
    if($a<$b){
	return $b;
    } else {
	return $a;
    }
}
