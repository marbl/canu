#!/usr/local/bin/perl

my %i2u = ();

$readIID2UIDmap = shift;
if($readIID2UIDmap eq ""){
    print STDERR "This script requires read IID-to-UID map argument!\n";
    exit(-1);
}

$cutoff = shift;

if($cutoff eq ""){
  $cutoff = 10000;
}

open(I2U,$readIID2UIDmap);
while(<I2U>){
    @w=split;
    $i2u{$w[0]}=$w[1];
}
close(I2U);

while(<>){
  chomp;
  $tag =substr($_,0,4);
  if($tag eq "{UTG"){
    $inIUM=1;
    next;
  }
  if($inIUM){
    if($tag eq "acc:"){
      $id = substr($_,5,length($_)-6);
      ($uid,$readiid) = split /,/,$id;
      next;
    }
    if($tag eq "cns:"){
      $inSeq=1;
      $seq="";
      next;
    }
    if($inSeq){
      if($_ eq "."){
	if(length($seq) >= $cutoff){
	  print ">$uid /readIID=$readiid /readUID=$i2u{$readiid}\n$seq\n";
	}
	$inSeq=0;
	$inIUM=0;
      } else {
	tr /-//d;
	$seq .= $_;
      }
    }
  }
}


