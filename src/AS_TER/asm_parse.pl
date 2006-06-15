#!/usr/local/bin/perl

use FileHandle;
#use FileUtil;


$frgplace = shift;
if($frgplace eq ""){
  print STDERR "ERROR: We require at least a prefix of file names\n";
  exit(-1);
}

$ctglens = shift;

if($ctglens eq ""){
  $prefix=$frgplace;
  $frgplace=$prefix . ".frgctg";
  $ctglens=$prefix . ".ctglen";
  $ctgscf=$prefix . ".ctgscf";
  $scflens=$prefix . ".scflen";
  $inputname = $prefix . ".asm";
} else {
  $ctgscf = shift;
  $scflens = shift;
  $inputname = shift;
}

$fh = new FileHandle();
if($inputname eq ""){
  $fh = STDIN;
} else {
  $fh->open($inputname);
}

if($scflens eq ""){
  print STDERR "USAGE: asm_parse.pl frgplace ctglens ctgscf scflen\n";
  print STDERR "where\n";
  print STDERR "\tfrgplace = file to contain fragment on contig info\n";
  print STDERR "\tctglens = file to contain contig lengths\n";
  print STDERR "\tctgscf = file to contain contig on scaffold info\n";
  print STDERR "\tscflens = file to contain scafafold lengths\n";
  exit(-1);
}

open(FRGPLACE,">$frgplace");
open(CTGLENS,">$ctglens");
open(CTGSCF,">$ctgscf");
open(SCFLENS,">$scflens");

$minAsmGap = 50;

while(<$fh>){
  chomp;
  $tag = substr($_,0,4);
  if($tag eq "{CCO"){
    $incontig=1;
    $inscf=0;
    $indsc=0;
    next;
  }
  if($tag eq "{SCF"){
    $incontig=0;
    $inscf=1;
    $indsc=0;
    next;
  }
  if($tag eq "{SLK"){
      $inscf=0;
      $incontig=0;
      $indsc=0;
      next;
  }
  if($tag eq "{DSC"){
    $incontig=0;
    $inscf=0;
    $indsc=1;
    next;
  }
  if($incontig){
    if($tag eq "acc:"){
      $postfix = substr($_,4);
      $postfix =~ /\((\d+),(\d+)\)/;
      $ctgID= $1;
    }
    if ( $tag eq "cns:") {
      $inseq=1;
      $ctgLen=0;
      undef $ctgCoords;
      $ctgCoords[0]=0;
      $I = 0;
      next;
    }
    if ($inseq) {
      if ($_ eq ".") {
	$inseq=0;
	if ($ctgLen>0) {
	  print CTGLENS "$ctgID $ctgLen\n";
	  $contigLength{$ctgID} = $ctgLen;
	}
	next;
      } else {
	@c = split //,$_;
	$n = @c;
	#	print "Split $_ into $n chars\n";
	for ($i=0;$i<$n;$i++) {
	  $I++;
	  if ($c[$i] ne "-") {
	    $ctgLen++;
	  }
	  $ctgCoords[$I] = $ctgLen;
	}
      }
      next;
    }
    if ($tag eq "{MPS") {
      $inmps=1;
      next;
    }
    if ($inmps){
      if($tag eq "mid:") {
	$frgID = substr($_,4);
	next;
      }
      if ($tag eq "pos:") {
	$coords = substr($_,4);
	$coords =~ /(\d+),(\d+)/;
	#	print "parsed $coords into $1 and $2\n";
	$b = $ctgCoords[$1];
	$e = $ctgCoords[$2];
	if($1>$2){ # do not use $b and $e -- it is possible, albeit
	           # very unusual, for a fragment to fit entirely
	           # in a gap in the consensus seq, in which case
	           # b == e; we need $rev to be set right in this case ...
	  $t = $b;
	  $b = $e;
	  $e = $t;
	  $rev = 1;
	} else {
	  $rev = 0;
	}
	print FRGPLACE "$frgID $ctgID $b $e $rev\n";
	$inmps=0;
	next;
      }
    }
  }
  if($inscf){
    if($tag eq "acc:"){
      if($scfLen > 0){
	print SCFLENS "$scfID $scfLen\n";
      }
      $postfix = substr($_,4);
      $postfix =~ /\((\d+),(\d+)\)/;
      $scfID= $1;
      $scfLen = 0;
      next;
    }
    if($tag eq "noc:"){
      $noc = substr($_,4);
      $cpn = 0;
      next;
    }
    if($tag eq "{CTP"){
      $cpn++;
      next;
    }
    if ($tag eq "ct1:") {
      $ct1id = substr($_,4);
      next;
    }
    if ($tag eq "ct2:") {
      $ct2id = substr($_,4);
      next;
    }
    if ($tag eq "mea:") {
      $mean = int(substr($_,4));
      if ($mean < $minAsmGap) {
	$mean = $minAsmGap;
      }
      next;
    }
    if ($tag eq "ori:") {
      $ori = substr($_,4);
      if ($cpn == 1) {
	$ctgB = $scfLen;
	$scfLen +=  $contigLength{$ct1id};
	$ctgE = $scfLen;
	if ($ori eq "A" || $ori eq "O") {
	  $rev = 1;
	} else {
	  $rev = 0;
	}
	print CTGSCF "$ct1id $scfID $ctgB $ctgE $rev\n";
      }
      if ($cpn>1|| $noc>0) {
	$scfLen += $mean;
	$ctgB = $scfLen;
	$scfLen += $contigLength{$ct2id};
	$ctgE = $scfLen;
	if ($ori eq "A" || $ori eq "I") {
	  $rev = 1;
	} else {
	  $rev = 0;
	}
	print CTGSCF "$ct2id $scfID $ctgB $ctgE $rev\n";
      }
      next;
    }
  }
  if($indsc){
      if($tag eq "acc:"){
	  $dscID = substr($_,4);
	  next;
      }
      if($tag eq "ctg:"){
	  $ctgID = substr($_,4);
	  # since processScaffolds uses the contig UID for the dregs file ...
	  print CTGSCF "$ctgID $ctgID 0 $contigLength{$ctgID} 0\n";
#	  print CTGSCF "$ctgID $dscID 0 $contigLength{$ctgID} 0\n";

	  print SCFLENS "$ctgID $contigLength{$ctgID}\n";
	  $indsc=0;
	  next;
      }
  }
}
  
if($scfLen > 0){
  print SCFLENS "$scfID $scfLen\n";
}
