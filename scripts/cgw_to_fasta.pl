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

# scan commandline
use Cwd;
$basename = cwd();
use Getopt::Std;
(getopts('f:') and $opt_f) or die "Wrong argument list : call cgw_to_fasta.pl -f <path to inputs>\n";

$basename = $opt_f;

print "basename = $basename\n";

$cat_cmd = "cat $basename.cgw $basename.cgw_contigs $basename.cgw_scaffolds > $basename.cgw_total";
print "Executing: $cat_cmd\n";
system($cat_cmd);
$consensus_cmd = "consensus -P $basename.frgStore $basename.cgw_total";
print "Executing: $consensus_cmd\n";
system($consensus_cmd);
$ps_cmd = "process_scaffolds -Z -F $basename.frgStore -f $basename.fasta < $basename.cns ";
print "Executing: $cps_cmd\n";
system($ps_cmd);
print "All Done!\n";

