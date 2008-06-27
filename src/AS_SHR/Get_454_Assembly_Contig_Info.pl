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

# $Id: Get_454_Assembly_Contig_Info.pl,v 1.2 2008-06-27 06:29:21 brianwalenz Exp $

###############################################################################

use strict;
use Getopt::Std;

my $usage = "usage:
$0 <fasta files...>

	This program will:
		Extract data in the fasta file from a 454 assembly.
";


my $opt_f;

if($#ARGV==-1){
	die $usage;
}

while($opt_f=shift){
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

}
print STDERR "Completed.\n";

###############################################################################

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	#.........................................................................
	my $length=grab_from_defline($defline, "length");
	my $num_reads=grab_from_defline($defline, "numReads");
	if(!defined($num_reads)){
		$num_reads=1;
	}
	my $id;
	if($defline=~/^>(\S+)/){
		$id=$1;
	}

	#.........................................................................
	# Output sequence after it has been formatted to the specified width
	print STDOUT "$id\t$length\t$num_reads\n";
}

#------------------------------------------------------------------------------

sub grab_from_defline{
	my $defline=shift;
	my $key=shift;
	my $value;

	$defline=~/$key=(\S+)/;
	$value=$1;
	return($value);
}

#------------------------------------------------------------------------------

