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
#        Script:  concat-uoms.pl
#   Description:  Perl script that concatenates the results of
#                 extract-lsf-uoms.pl into a cgb file
#                 
#    Programmer:  Knut Reinert
#       Written:  30 Nov 99
#********************************************************************

$file = shift(@ARGV);
$blocksize  = shift(@ARGV);
$upperBound = shift(@ARGV);


print "concat-uoms.pl with parameters\n $file\n $blocksize\n $upperBound\n";


$i = 0;
system "rm $file.filtered.cgb";
system "rm stats/$file.O.cgm";
system "rm stats/$file.M.cgm";
system "rm stats/$file.I.cgm";
system "rm stats/$file.Y.cgm";
system "rm stats/$file.X.cgm";
system "rm stats/$file.R.cgm";


while( $i*$blocksize < $upperBound )
  {
    system "cat uoms/$file.$i >> $file.filtered.cgb";
    system "cat uoms/$file.O.stat.$i >> stats/$file.O.cgm";
    system "cat uoms/$file.M.stat.$i >> stats/$file.M.cgm";
    system "cat uoms/$file.I.stat.$i >> stats/$file.I.cgm";
    system "cat uoms/$file.Y.stat.$i >> stats/$file.Y.cgm";
    system "cat uoms/$file.X.stat.$i >> stats/$file.X.cgm";
    system "cat uoms/$file.R.stat.$i >> stats/$file.R.cgm";

    $i++;
  }
