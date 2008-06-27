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

# $Id: Filter_Short_Reads.pl,v 1.2 2008-06-27 06:29:21 brianwalenz Exp $

###############################################################################

use strict;
use Getopt::Std;
use vars qw($opt_f $opt_l);

my $MIN_LENGTH=64;

getopts("f:l:");
my $usage = "usage:
$0
	-f <fasta file>
	-l <length, optional>

	This program will:
		Will read in a FASTA file, and filter out reads <= $MIN_LENGTH.
";

if(!(
	defined($opt_f))){
	die $usage;
}

if(defined($opt_l)){
	$MIN_LENGTH=$opt_l;
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

sub process_record{
	my $defline = shift;
	my $sequence = shift;

	#.........................................................................
	# Output sequence after it has been formatted to the specified width
	my $length=length($sequence);
	if($length>$MIN_LENGTH){
	    print STDOUT "$defline\n";
	    my $width=50;
	    my $pos=0;
	    do{
		    my $out_width=($width>$length)?$length:$width;
		    print STDOUT substr($sequence, $pos, $width) . "\n";
		    $pos+=$width;
		    $length-=$width;
	    }while($length>0);
	}
}

#------------------------------------------------------------------------------

