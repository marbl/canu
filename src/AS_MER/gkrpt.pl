#!/usr/bin/perl

use strict;

my %kmers;

while (!eof(STDIN)){
  my $count = <STDIN>;  chomp $count;
  my $kmer  = <STDIN>;  chomp $kmer;
  $count =~ s/>//g;
  $count = int($count);
  $kmer =~ tr/acgt/ACGT/;
  my $numbad = $kmer =~ tr/ACGT//c;
  if ($numbad > 0){
    die "$numbad nonacgtACGT characters in kmer $kmer\n";
  }
  if (defined($kmers{$kmer})){
    die "kmer repeated in input $kmer\n";
  }
  $kmers{$kmer} = $count;
}

foreach my $kmer (sort { $kmers{$b} <=> $kmers{$a} } (keys %kmers)){
  my $startkmer = my $curkmer = $kmer;
  my $startcount = my $curcount = $kmers{$kmer};
  if ($curcount < 0) {
    next;
  } else {
    $kmers{$kmer} = - $curcount;
  }
  my $currpt = $curkmer;
  for ( ; ; ) {
    my $nextcount = -1;
    my $realcount = -1;
    my $maxkmer;
    my $realkmer;
    my $tmpkmer;
    my $maxcount = -1;
    $curkmer = (substr $curkmer, 1) . "A";
    my $rckmer = reverse $curkmer;
    $rckmer =~ tr/ACGT/TGCA/;
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    substr($curkmer, -1) = "C";
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      substr($rckmer, 0, 1) = "G";
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    substr($curkmer, -1) = "G";
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      substr($rckmer, 0, 1) = "C";
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    substr($curkmer, -1) = "T";
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      substr($rckmer, 0, 1) = "A";
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    if (($realcount < 0) || ($realcount < ($curcount / 2))) {
      last;
    } else {
      $curkmer = $maxkmer;
      $curcount = $realcount;
      $kmers{$realkmer} = - $realcount;
      $currpt .= (substr $curkmer, -1);
    }
  }
  $curcount = $startcount;
  $curkmer = $startkmer;

  for ( ; ; ) {
    my $nextcount = -1;
    my $realcount = -1;
    my $maxkmer;
    my $realkmer;
    my $tmpkmer;
    my $maxcount = -1;
    $curkmer = "A" . (substr $curkmer, 0, -1);
    my $rckmer = reverse $curkmer;
    $rckmer =~ tr/ACGT/TGCA/;
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    substr($curkmer, 0, 1) = "C";
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      substr($rckmer, -1) = "G";
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    substr($curkmer, 0, 1) = "G";
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      substr($rckmer, -1) = "C";
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    substr($curkmer, 0, 1) = "T";
    if (defined($kmers{$curkmer})){
      $nextcount = $kmers{$curkmer};
      $tmpkmer = $curkmer;
    } else {
      substr($rckmer, -1) = "A";
      if (defined($kmers{$rckmer})){
	$nextcount = $kmers{$rckmer};
	$tmpkmer = $rckmer;
      } else {
	$nextcount = -1;
      }
    }
    if ((abs $nextcount) > $maxcount){
      $maxcount = abs $nextcount;
      $realcount = $nextcount;
      $realkmer = $tmpkmer;
      $maxkmer = $curkmer;
    }
    if (($realcount < 0) || ($realcount < ($curcount / 2))) {
      last;
    } else {
      $curkmer = $maxkmer;
      $curcount = $realcount;
      $kmers{$realkmer} = - $realcount;
      $currpt = (substr $curkmer, 0, 1) . $currpt;
    }
  }
  if ((my $lenrpt = length $currpt) > $ARGV[0]) {
    print ">$startkmer $startcount $lenrpt\n$currpt\n";
  }
}
