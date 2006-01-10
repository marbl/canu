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

# $Id: Report_Duplicate_Statistics.pl,v 1.1 2006-01-10 22:44:42 kli1000 Exp $

###############################################################################

use strict;
use Getopt::Std;

print "Started...\n";

my $usage = "usage: 
$0 <fasta files...>
	This program will:
		Will read in a FASTA file, and print out id's that have the same sequence.
";

if($#ARGV==-1){
	die "$usage\n";
}

my %seq_hash;
my $num_reads=0;

my $fname;
while($fname=shift){

	print "\n\n$fname\n";
	open(FASTA_FH, "<$fname") || die "Could not open $fname\n";

	###############################################################################
	%seq_hash=();
	$num_reads=0;

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

	my $dup_counter=0;
	my @dupl_histogram;
	my $i;
	my $max_dups=0;

	for($i=2; $i<1000; $i++){
		$dupl_histogram[$i]=0;
	}

	print "\nNum Reads\t$num_reads\n\n";

	foreach my $seq(keys %seq_hash){

		my $num_dups=$#{$seq_hash{$seq}}+1;
		if($num_dups>1){
			$dupl_histogram[$num_dups]++;
			if($max_dups<$num_dups){
				$max_dups=$num_dups;
			}
		}
	}

	my $num_groups=0;
	for($i=2; $i<=$max_dups; $i++){
		my $normalized=$dupl_histogram[$i]/$num_reads;
		printf "$i\t$dupl_histogram[$i]\t%3.6f\n", $normalized;
		$num_groups+=$dupl_histogram[$i];
	}

	printf "Num Groups:\t$num_groups\t%3.6f\n", $num_groups/$num_reads;


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

	push @{$seq_hash{$sequence}}, $id;
	$num_reads++;

}

#------------------------------------------------------------------------------

