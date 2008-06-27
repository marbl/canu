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

# $Id: Extract_Reads_From_ACE.pl,v 1.2 2008-06-27 06:29:21 brianwalenz Exp $

use strict;
use Getopt::Std;
use FileHandle;
use vars qw($opt_a $opt_f $opt_s);

my $MIN_COVERAGE=2;
my $MIN_CONTIG_SIZE=600;

getopts("a:f:s:");
my $usage = "usage:
$0
        -a <ace file, input>
        -f <fasta file, output>

	This program will extract reads from an ace file.

";

if(!(
        defined($opt_a) &&
        defined($opt_f))){
        die $usage;
}

###############################################################################

my $fh=new FileHandle "<$opt_a";
my $fasta_fh=new FileHandle ">$opt_f";

# || die "Could not open the 454 ace file: $opt_a\n";


my ($num_contigs, $num_reads)=read_AS($fh);

print STDERR "Number of Contigs: $num_contigs\n";
print STDERR "Number of Reads:   $num_reads\n";

my $contig_idx;
for($contig_idx=0; $contig_idx<$num_contigs; $contig_idx++){

	my %read_length_hash;
	my %read_position_hash;

	my ($contig_id, $num_consensus_bases, $num_reads, $num_segments, $complementation, $consensus_sequence)=
		read_CO($fh);

	#print STDOUT "$contig_id: num_reads: $num_reads \n$consensus_sequence\n";

	my $quality=read_BQ($fh);

	my $read_idx;
	for($read_idx=0; $read_idx<$num_reads; $read_idx++){
		my ($read_id, $complementation, $consensus_start_pos)=read_AF($fh);
		$read_position_hash{$read_id}=$consensus_start_pos-1; #convert to space based
	}

	my ($base_line_start, $base_line_end, $base_line_read_id)=read_BS($fh);

	for($read_idx=0; $read_idx<$num_reads; $read_idx++){
		my ($read_id, $num_padded_bases, $num_read_info_items, $num_read_tags, $read_sequence)=
			read_RD($fh);

		if(!($read_id=~/fake/)){
			$read_sequence=~s/\*//g;
			print $fasta_fh ">$read_id\n$read_sequence\n";
		}

		my ($qual_start, $qual_end, $align_start, $align_end)=read_QA($fh);
		my ($null)=read_DS($fh);
	}


}

###############################################################################
###############################################################################

sub read_AS{
	my $fh=shift;
	my ($id, $num_contigs, $num_reads);

	while(<$fh>){
		chomp;
		($id, $num_contigs, $num_reads)=split /\s+/;
		if($id eq "AS"){
			return ($num_contigs, $num_reads);
		}
	}
	die "Could not find AS to read.\n";
}

###############################################################################

sub read_CO{
	my $fh=shift;

	my ($id, $contig_id, $num_bases, $num_reads, $num_segments, $complementation, $sequence);

	while(<$fh>){
		chomp;
		($id, $contig_id, $num_bases, $num_reads, $num_segments, $complementation, $sequence)=
			split /\s+/;

		if($id eq "CO"){
			while(<$fh>){
				chomp;
				if($_ eq ""){
				    last;
				}else{
				    $sequence.=$_;
				}
			}
			return($contig_id, $num_bases, $num_reads, $num_segments, $complementation, $sequence);
		}
	}
	die "Could not find CO to read.\n";
}

###############################################################################

sub read_BQ{
	my $fh=shift;

	my ($id, $sequence);

	while(<$fh>){
		chomp;
		($id)=split /\s+/;

		if($id eq "BQ"){
			while(<$fh>){
				chomp;
				if($_ eq ""){
					last;
				}else{
					$sequence.=$_;
				}
			}
			return($sequence);
		}
	}
	die "Could not find BQ to read.\n";
}

###############################################################################

sub read_AF{
	my $fh=shift;

	my ($id, $read_id, $complementation, $start);

	while(<$fh>){
		chomp;
		($id, $read_id, $complementation, $start)=split /\s+/;
		if($id eq "AF"){
			return($read_id, $complementation, $start);
		}
	}
	die "Could not find AF to read.\n";
}

###############################################################################

sub read_BS{
	my $fh=shift;

	my ($id, $start, $end, $read_id);

	while(<$fh>){
		chomp;
		($id, $start, $end, $read_id)=split /\s+/;
		if($id eq "BS"){
			return($start, $end, $read_id);
		}
	}
	die "Could not find BS to read.\n";
}

###############################################################################

sub read_RD{
	my $fh=shift;

	my ($id, $read_id, $num_bases, $num_read_info_items, $num_read_tags, $sequence);

	while(<$fh>){
		chomp;
		my ($id, $read_id, $num_bases, $num_read_info_items, $num_read_tags)=
			split /\s+/;
		if($id eq "RD"){
			while(<$fh>){
				chomp;
				if($_ eq ""){
					last;
				}else{
					$sequence.=$_;
				}
			}
			return($read_id, $num_bases, $num_read_info_items, $num_read_tags, $sequence);
		}
	}
	die "Could not find RD to read.\n";
}

###############################################################################

sub read_QA{
	my $fh=shift;

	my ($id, $qual_start, $qual_end, $clip_start, $clip_end);

	while(<$fh>){
		chomp;
		my ($id, $qual_start, $qual_end, $clip_start, $clip_end)=split /\s+/;
		if($id eq "QA"){
			return($qual_start, $qual_end, $clip_start, $clip_end);
		}
	}
	die "Could not find QA to read.\n";
}

###############################################################################

sub read_DS{
	my $fh=shift;
	my $id;
	while(<$fh>){
		chomp;
		my ($id)=split /\s+/;
		if($id eq "DS"){
			return("not implemented");
		}
	}
	die "Could not find DS to read.\n";
}











