#!/usr/bin/perl 
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
###########################################################################
#
#********************************************************************
#        Script:  extract-lsf-uoms.pl
#   Description:  Perl script that invokes calls of extract-quality-uoms
#                 submits the jobs to the lsf queue assembly-filter
#                 
#    Programmer:  Knut Reinert
#       Written:  30 Nov 99
#********************************************************************

$queue = shift(@ARGV);
# name of the lsf queue
$file  = shift(@ARGV);
# file prefix 
$blocksize  = shift(@ARGV);
# how many pieces each job (rec. for DROS 2000000)
$upperBound = shift(@ARGV);
# upper bound for the number of UOM messages
$threshold  = shift(@ARGV);
# filtering threshild (DROS runs 0.5)
$allocate   = shift(@ARGV);
# set to the number of IUMs

if( $allocate == 0 )
  {
    $allocate = 1024;
  }

print "\ncalling extract_uoms.pl with parameters\n $queue\n $file\n $blocksize\n $upperBound\n $threshold\n $allocate\n";

system "mkdir uoms";

$i = 0;

while( $i*$blocksize < $upperBound )
  {
    if( ! -e "uoms/$file.$i" )
      {
	if( $i == 0 )
	  {
	    system "bsub -C 0 -R \"select[mem>4000]rusage[mem=4000:swap=5000:duration=30]\" -q $queue -o uoms/$file.$i.log extract-quality-uoms -r -b $blocksize -n $i  -t $threshold -a $allocate -f $file &";
	  }
	else
	  {
	    system "bsub -C 0 -R \"select[mem>4000]rusage[mem=4000:swap=5000:duration=30]\" -q $queue -o uoms/$file.$i.log extract-quality-uoms -b $blocksize -n $i -t $threshold -a $allocate -f $file &";
	  }
      }
    $i++;
  }

