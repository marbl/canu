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

$file  = shift(@ARGV);
$blocksize  = shift(@ARGV);
$upperBound = shift(@ARGV);
$threshold  = shift(@ARGV);
$allocate   = shift(@ARGV);

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
	    system "extract-quality-uoms -r -b $blocksize -n $i  -t $threshold -a $allocate -f $file &";
	  }
	else
	  {
	    system "extract-quality-uoms -b $blocksize -n $i -t $threshold -a $allocate -f $file &";
	  }
      }
    $i++;
  }



