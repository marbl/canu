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
#        Script:  overlap_script.pl
#   Description:  Perl script to save a batch of new fragments in
#                 the fragment store and then compute all overlaps
#                 between *all* fragments in the store and the new
#                 fragments
#                 
#    Programmer:  Art. Delcher
#   Stolen from:  Knut Reinert
#       Revised:  11 April 2000
#********************************************************************

$lsfScriptName = "fragovl.script";

# command names and error messages 
%cmd = (
'populator'       => '$AS_BIN/PopulateFragStore',
'create-lsf-jobs' => '$AS_BIN/make-ovl-script',
'run-lsf-jobs'    => './' . "$lsfScriptName"
);


%error = (
'populator'       => "populator failed. please check log file and $AS_BIN!\n",
'create-lsf-jobs' => "create-lsf-jobs failed. please check log file !\n",
'run-lsf-jobs'    => "run-lsf-jobs failed. please check log file !\n",
'rename'          => "rename filed, check log file !\n"
);


# start script
$start = (times)[3];

use Cwd;
$cwd = cwd();

# scan commandline
use Getopt::Std;
getopts('1aAcefhi:K:l:o:Pp:q:Qt:vwX')
    or die "Wrong argument list : call overlap_script.pl [-1aAcfhPpQvwX]\
  [-p <parameter file>] \
  [-i <inputStore>] [-K <kmerfreqlimit>] [-o <outputStore>]\
  [-l <olaplimit>] [-q <lsfqueue>]\
  [-t <numthreads>] <inputFile>.<ext>\n";

if( $opt_h == 1 )
{
    print "call overlap_script.pl [-1aAcfhpPQvw] [-i <inputStore>] [-o <outputStore>]\n";
    print "   [-l <olaplimit>] [-q <lsfqueue>] [-t <numthreads>] <inputFile>.<ext>\n";
    print "Updates fragment store and then computes overlaps between\n";
    print "  old fragments and those newly added\n";
    print "Outputs:\n";
    print "  <inputFile>.ovl    filtered input messages and overlap messages\n";
    print "  <inputFile>.range  range of new fragment numbers added to store\n";
    print "Options:\n";
    print "  -1  create a single .ovl file (instead of separate ones\n";
    print "      for each LSF job\n";
    print "  -a  append with backup (requires -i)\n";
    print "  -A  append with no backup (requires -i)\n";
    print "  -c  create (requires -o)\n";
    print "  -e  use exact-multiple batch sizes for LSF overlap jobs\n";
    print "  -f  force.  If output store exists, nuke it.\n";
    print "  -h  print this message\n";
    print "  -i  specifies input fragStore\n";
    print "  -K  specifies min hash-table kmer frequency to ignore\n";
    print "  -l  specifies max olaps per old frag end per hash batch\n";
    print "  -o  specifies output fragStore\n";
    print "  -P  proto (ASCII) output (default is binary)\n";
    print "  -p  specifies the parameters file.\n";
    print "  -q  specifies LSF queue for jobs\n";
    print "  -Q  do NOT use LSF for overlap jobs\n";
    print "  -t  specifies number of parallel threads for overlap jobs\n";
    print "  -v  verbose mode\n";
    print "  -w  filter overlaps with too many errors in a window\n";
    print "  -X  activate the expert options\n";
    exit;
}

if(@ARGV[0] eq "")
{
    print "No input file specified\n";
    exit;
}

$inputFileName = @ARGV[0];

$pos = -1;
while (($place = index($inputFileName, ".", $pos + 1)) > -1)
{
    $pos = $place;
}

if($pos == -1)
{
    $prefix = $inputFileName;
}
else
{
    $prefix = substr($inputFileName, 0, $pos);
}

print "\n*** Running overlap_script\n";


# Run the populator

$populatorOptions = "";
if( $opt_a == 1)
{
    $populatorOptions = $populatorOptions . " -a";
}
if( $opt_A == 1)
{
    $populatorOptions = $populatorOptions . " -A";
}
if( $opt_c == 1)
{
    $populatorOptions = $populatorOptions . " -c";
}
if( $opt_f == 1)
{
    $populatorOptions = $populatorOptions . " -f";
}
if( $opt_i ne "")
{
    $populatorOptions = $populatorOptions . " -i $opt_i";
}
if( $opt_o ne "")
{
    $populatorOptions = $populatorOptions . " -o $opt_o";
}
if( $opt_P == 1)
{
    $populatorOptions = $populatorOptions . " -P";
}
if( $opt_p ne "")
{
    print "\n*** -p option is out of commission.\n";
}
if( $opt_v == 1)
{
    $populatorOptions = $populatorOptions . " -v";
}
if( $opt_X == 1)
{
    $populatorOptions = $populatorOptions . " -X";
}

print "populator options = $populatorOptions\n";
print "output to file $prefix.ovlTMP\n";

if( 0 != system "$cmd{'populator'} $populatorOptions -V $prefix.ovlTMP $inputFileName")
  { die $error{'populator'}; }

if(0 != system "mv $prefix.ovlTMP $prefix.ovl")
{
    die $error{'rename'};
}


# Build the overlapper lsf script

if(  $opt_o ne "")
{
    $fragStore = $opt_o;
}
else
{
    $fragStore = $opt_i;
}

$overlapperOptions = "";
if( $opt_1 == 1)
{
    $overlapperOptions = $overlapperOptions . " -1";
}
if( $opt_e == 1)
{
    $overlapperOptions = $overlapperOptions . " -e";
}
if( $opt_K ne "")
{
    $overlapperOptions = $overlapperOptions . " -K $opt_K";
}
if( $opt_l ne "")
{
    $overlapperOptions = $overlapperOptions . " -l $opt_l";
}
if( $opt_P == 1)
{
    $overlapperOptions = $overlapperOptions . " -P";
}
if( $opt_q ne "")
{
    $overlapperOptions = $overlapperOptions . " -q $opt_q";
}
if( $opt_Q == 1)
{
    $overlapperOptions = $overlapperOptions . " -Q";
}
if( $opt_t ne "")
{
    $overlapperOptions = $overlapperOptions . " -t $opt_t";
}
if( $opt_w == 1)
{
    $overlapperOptions = $overlapperOptions . " -w";
}

$overlapperOptions = $overlapperOptions . " -o $prefix -s $lsfScriptName";

print "overlapper options = $overlapperOptions\n";


if( 0 != system "$cmd{'create-lsf-jobs'} $overlapperOptions $fragStore $prefix.range")
  { die $error{'create-lsf-jobs'}; }


# Run the overlapper lsf script

if( 0 != system "$cmd{'run-lsf-jobs'}")
  { die $error{'run-lsf-jobs'}; }


# Rename the .ovl output file back to its real name

#if  0
if(0 != system "mv $prefix.ovlTMP $prefix.ovl")
{
    die $error{'rename'};
}
#endif

# if we are not dead here we succesfully processed the batch
print "\n############# DONE ################################ \n";
print "+++ successfully processed file $inputFileName +++\n";
$end = (times)[3];
printf "processing took  %.2f system seconds\n",$end -$start;
