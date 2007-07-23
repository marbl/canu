#!/usr/local/bin/perl

##########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received (LICENSE.txt) a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
##########################################################################

# $Id: FASTA_to_frg_file.pl,v 1.5 2007-07-23 17:23:41 brianwalenz Exp $

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_q $opt_f);

print STDERR "WARNING!\n";
print STDERR "WARNING!\n";
print STDERR "WARNING!  This is UNTESTED frag format version 2 code!  It will work\n";
print STDERR "WARNING!  ONLY with assemblers more recent than May 10, 2007.  (If it\n";
print STDERR "WARNING!  tests ok, kindly remove this warning, thanks.)\n";
print STDERR "WARNING!  Complaints to Bri.\n";
print STDERR "WARNING!\n";
print STDERR "WARNING!\n";

my $DEFAULT_QUAL=20;
my $LOW_QUAL_DIVISOR=4;

getopts("f:q:");
my $usage = "usage:
$0
        -f <fasta file>
	-q <quality values, default $DEFAULT_QUAL>

	This program will read in a multi-FASTA file and generate
	a Celera Assembler FRG file.

";

if(!(
        defined($opt_f))){
        die $usage;
}

if(defined($opt_q)){
	$DEFAULT_QUAL=$opt_q;
}
print STDERR "Default Quality Value: $DEFAULT_QUAL\n";

my $filename=$opt_f;
my $time=time;

###############################################################################
# Print Library Record

my $uidServ = new Annotation::UID(2, "seq454");
my $libId   = $uidServ->incrUID;

print STDOUT "{VER\n";
print STDOUT "ver:2\n";
print STDOUT "}\n";
print STDOUT "{LIB\n";
print STDOUT "act:A\n";
print STDOUT "acc:$libId\n";
print STDOUT "ori:U\n";
print STDOUT "mea:0.0\n";
print STDOUT "std:0.0\n";
print STDOUT "src:\n";
print STDOUT ".\n";
print STDOUT "nft:1\n";
print STDOUT "doNotOverlapTrim=1\n";
print STDOUT ".\n";
print STDOUT "}\n";

my $frag_id_counter=1;

###############################################################################

my $fasta_fh=new FileHandle "<$filename";

print STDERR "Processing FASTA file...\n";
my ($defline, $prev_defline, $sequence);
while(<$fasta_fh>){
	chomp;
	
	if(/^>/){
		$defline=$_;
		if($sequence ne ""){
			process_record($prev_defline, $sequence);
			$sequence="";
		}
		$prev_defline=$defline;
	}else{
		$sequence.=$_;
	}
}
process_record($prev_defline, $sequence);

print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $starting_length=length($sequence);

	####################################################################

	# Generate quality values
	my $qual_vals=chr($DEFAULT_QUAL + ord("0")) x $starting_length;

	# Reduce quality values of lowercase nucleotides
	my @nucs=split //, $sequence;
	my @quals=split //, $qual_vals;
	my $i;
	for($i=0; $i<=$#nucs; $i++){
		if($nucs[$i] eq lc($nucs[$i])){
			$quals[$i]=chr($DEFAULT_QUAL/$LOW_QUAL_DIVISOR + ord("0"));
		}	
	}

	$qual_vals=join "", @quals;
	$sequence=join "", @nucs;

	####################################################################

	my $frag_id;
	if($defline=~/^>(\S+)/){
		$frag_id=$1;
	}

	print STDOUT "{FRG\n";
	print STDOUT "act:A\n";
	print STDOUT "acc:$frag_id\n";
	print STDOUT "rnd:1\n";
	print STDOUT "sta:G\n";
	print STDOUT "lib:$libId\n";
	print STDOUT "pla:0\n";
	print STDOUT "loc:0\n";
	print STDOUT "src:\n.\n";

	print STDOUT "seq:\n";
	
	# Output sequence
	my $length=$starting_length;
	my $width=70;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT lc(substr($sequence, $pos, $width)) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

	print STDOUT ".\n";
	print STDOUT "qlt:\n";

	# Output quality values
	$length=$starting_length;
	$pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr($qual_vals, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

        print STDOUT "hps:\n.\n";

	print STDOUT ".\n";
	print STDOUT "clr:0,$starting_length\n";
	print STDOUT "}\n";

	$frag_id_counter++;
}

#------------------------------------------------------------------------------
