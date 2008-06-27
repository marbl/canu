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

# $Id: Shred_Contigs.pl,v 1.5 2008-06-27 06:29:21 brianwalenz Exp $

use strict;
use Getopt::Std;
use FileHandle;
use FindBin qw($Bin);
use lib $Bin;
use Annotation::UID;
use vars qw($opt_r $opt_f);

getopts("f:r:");
my $usage = "usage:
$0
        -f <fasta file>
	-r <read length, default 600>

	This program will read in a multi-FASTA file and generate shredded
	simulated reads for each record depending on the depth of coverage
	specified in the defline.

	The tag /average_coverage, should exists already in the defline.

";

if(!(
        defined($opt_f))){
        die $usage;
}

my $read_length=$opt_r;
if(!defined($read_length)){
	$read_length=600;
}
print STDERR "Read Length: $read_length\n";

my $logConf = q(
log4perl.category.GUID          = WARN, Screen
log4perl.appender.Screen        = Log::Log4perl::Appender::Screen
log4perl.appender.Screen.stderr = 0
log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
);
Log::Log4perl->init(\$logConf);
my $uidBatchSize = 1000;
my $uidNamespace = 'seq454';
my $uidServ = new Annotation::UID( $uidBatchSize, $uidNamespace);

###############################################################################

my $fasta_fh=new FileHandle "<$opt_f";

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

	my $avg_cov=grab_from_defline($defline, "average_coverage");

	my $seq_id;
	if($defline=~/^>(\S+)/){
		$seq_id=$1;
	}

	my $seq_len=length($sequence);
	my ($begin_shred_ref, $end_shred_ref)=shred($seq_len, $avg_cov, $read_length);
	my $num_shreds=$#{$begin_shred_ref}+1;
	my $accomplished_coverage=$num_shreds*$read_length/$seq_len;

	# Output sequence after it has been formatted to the specified width
	my $shred_idx;
	for($shred_idx=0; $shred_idx<$num_shreds; $shred_idx++){

		my $shredded_sequence=substr(
			$sequence,
			${$begin_shred_ref}[$shred_idx],
			${$end_shred_ref}[$shred_idx]-${$begin_shred_ref}[$shred_idx]);

        my $frgId = $uidServ->incrUID;
		print STDOUT
            ">$frgId " ,
			"/contig=$seq_id\.$shred_idx " ,
			"/target_coverage=$avg_cov " ,
			"/accomplished_coverage=$accomplished_coverage " ,
			"/input_length=$seq_len " ,
			"/range=${$begin_shred_ref}[$shred_idx]-" ,
			       "${$end_shred_ref}[$shred_idx]\n";

		my $length=length($shredded_sequence);
		my $width=50;
		my $pos=0;
		do{
			my $out_width=($width>$length)?$length:$width;
			print STDOUT substr($shredded_sequence, $pos, $width) . "\n";
			$pos+=$width;
			$length-=$width;
		}while($length>0);
	}
}

#------------------------------------------------------------------------------

sub shred{
	my $seq_len=shift;
	my $target_coverage=shift;
	my $read_len=shift;

	my @begins=();
	my @ends=();

	#
	#                  |*******|
	#                  |###############|
	# |-------------------------------------------------|
	#  ----------------1----------------
	#          ----------------2----------------
	#                  ----------------3----------------
	#
	#	#### represents the distance between center of read 1 and read 3
	#            [$center_range_width]
	#       **** represents the distance between centers of consective reads
	#            [$center_increments]
	#

    my $shred_len = $read_len;
    $shred_len = $seq_len - 50 if $seq_len < $read_len;

    my $num_reads=int($seq_len*$target_coverage/$shred_len);
    my $center_range_width=$seq_len-$shred_len;
    if($num_reads==1){
        push @begins, 0;
        push @ends, $shred_len;
    }else{
        my $center_increments=$center_range_width/($num_reads-1);

# Cap the number of reads we will make so that we don't get
# redundant reads

        my $i;
        my ($prev_begin, $prev_end)=(-1,-1);
        for($i=0; $i<$num_reads; $i++){
            my $begin=$center_increments*$i;
            my $end=$begin+$shred_len;

            $begin=int($begin);
            $end=int($end);
#print "$begin-$end\n";

            if($begin!=$prev_begin || $end!=$prev_end){
                push @begins, $begin;
                push @ends, $end;
                $prev_begin=$begin;
                $prev_end=$end;
            }
        }
    }

    return(\@begins, \@ends);
}


#------------------------------------------------------------------------------

sub grab_from_defline{
    my $defline=shift;
    my $key=shift;
    my $value;

	$defline=~/\/$key=(\S+)/;
	$value=$1;
	return($value);
}

#------------------------------------------------------------------------------


