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

# $Id: Almost_Randomly_Substitute_Ns.pl,v 1.1 2006-01-10 22:44:42 kli1000 Exp $

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f);

getopts("f:");
my $usage = "usage: 
$0 
	-f <fasta file>

	This program will:
		Will read in a FASTA file, and remove N's by substituting them with
		A/T/C/G.  The nucleotides are selected randomly based on what nucleotides
		are not neighboring the N on both sides.  Substituted N's will be in
		lowercase.  If there is a stretch of N's, the internal N's will be
		randomly selected based on all 4 nucleotides, not just on avoiding
		immediate neighbors.
";

if(!(
	defined($opt_f))){
	die $usage;
}

###############################################################################
# Make sure files open before wasting any time doing anything

open(FASTA_FH, "<$opt_f") || die "Could not open $opt_f\n";

###############################################################################

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

print STDERR "Completed.\n";

###############################################################################

sub substitute_N{

	my $left=shift;
	my $right=shift;

	my %choices;
	$choices{"A"}=1;
	$choices{"T"}=1;
	$choices{"G"}=1;
	$choices{"C"}=1;

	delete $choices{$left};
	delete $choices{$right};

	my @remaining=keys %choices;

	my $random=rand() * ($#remaining+1);
	return($remaining[int($random)]);

}

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	my @arr=split //, $sequence;
	my $i;
	my $Ns_found=0;
	for($i=0; $i<=$#arr; $i++){
		if($arr[$i] eq "N"){
			$arr[$i]=lc(substitute_N($arr[$i-1], $arr[$i+1]));
		}
	}

	$sequence=join "", @arr;

	#.........................................................................
	# Output sequence after it has been formatted to the specified width
	my $length=length($sequence);
	print STDOUT "$defline\n";
	my $width=60;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print STDOUT substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);
}

#------------------------------------------------------------------------------

