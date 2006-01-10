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

# $Id: Remove_Duplicate_Sequences.pl,v 1.1 2006-01-10 22:44:42 kli1000 Exp $

###############################################################################

use strict;
use Getopt::Std;

my $usage = "usage: 
$0 <fasta files...>
	This program will:

		Read in a multi-FASTA file, and generate a non redundant multi-FASTA file 
		where the defline of the unique sequence is the defline of the first 
		occurrence of the sequence.
";

if($#ARGV==-1){
	die "$usage\n";
}

my %seq_hash;
my $num_reads=0;

my $fname;
while($fname=shift){

	print "Working on: $fname\n";
	open(FASTA_FH, "<$fname") || die "Could not open $fname\n";
	open(OUT_FH, ">$fname\.nr") || die "Could not open $$fname\.nr";

	###############################################################################
	%seq_hash=();

	print STDERR "Processing FASTA file...\n";
	my ($defline, $prev_defline, $sequence);
	while(<FASTA_FH>){
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

	close(FASTA_FH);
	close(OUT_FH);

	print STDERR "Completed.\n";
}

###############################################################################


sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}

	if(!defined($seq_hash{$sequence})){
		print OUT_FH "$defline\n";

		my $length=length($sequence);
		my $width=50;
		my $pos=0;
		do{
			my $out_width=($width>$length)?$length:$width;
			print OUT_FH substr($sequence, $pos, $width) . "\n";
			$pos+=$width;
			$length-=$width;
		}while($length>0);

		$seq_hash{$sequence}=$id;
	}

}

#------------------------------------------------------------------------------

