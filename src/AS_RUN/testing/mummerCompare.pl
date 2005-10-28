#!/usr/local/bin/perl
#
###########################################################################
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
# $Id: mummerCompare.pl,v 1.3 2005-10-28 19:54:41 catmandew Exp $
#

#  Script to analyze a show-coords output file using the
#  corresponding assembly.scaff files as a guide
#
#  To use, run:
#    nucmer ref.scaffolds.fasta  query.scaffolds.fasta
#
#  nucmer creates out.cluster and out.delta.
#
#  Then run:
#    show-coords -THcl -I 99 out.delta > showCoordsResults.txt
#
#  Then run this script:
#     showCoordsCompare -r ref.scaff -q query.scaff -c showCoordsResults.txt
#
#
#   Written by Ian Dew
#

use Carp;
use strict;
use FileHandle;
use Getopt::Long;

my $MY_VERSION = " Version 1.01 (Build " . (qw/$Revision: 1.3 $/ )[1]. ")";
my $MY_APPLICATION = "showCoordsCompare";

my $HELPTEXT = qq~
Compare two caqc-generated qc files

    showCoordsCompare  [options]  -r ref.scaff -q query.scaff -c show-coordsFile

    -r ref.scaff        reference sequence .scaff file

    -q query.scaff      query sequence .scaff file

    -c show-coordsFile  output from show-coords -THcl -I 99

    options:
      -h               Print help.
  
      -v <level>       Set verbosity to level.

$MY_VERSION

~;


my $INTRA_SCAFFOLD_GAP_LENGTH = 100;

######################################################################
# Parse the command line
######################################################################
my $helpRequested;
my $verboseLevel = 0;
my $rfile;
my $qfile;
my $cfile;

GetOptions("r=s" => \$rfile,
           "q=s" => \$qfile,
           "c=s" => \$cfile,
           "h|help" => \$helpRequested,
           "v|V|verbose:1" => \$verboseLevel
           ) or die $HELPTEXT;

if($helpRequested)
{
  print STDERR "Help requested:\n\n";
  print STDERR $HELPTEXT;
  exit 0;
}

if(! $rfile || ! -f $rfile ||
   ! $qfile || ! -f $qfile ||
   ! $cfile || ! -f $cfile)
{
  print STDERR "Please specify valid input files\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}


######################################################################
# Read in .scaff files
######################################################################
my %rdata = ReadScaffFile($rfile);
DumpScaffData(%rdata);
my %qdata = ReadScaffFile($qfile);




######################################################################
# Stream through show-coords output file & analyze
######################################################################

my $ifh = new FileHandle $cfile, "r"
  or die "Failed to open $cfile for reading";
while(<$ifh>)
{
  s/[\n\r\cZ]//g;

  my @fields = split;

  # coords are 1-offset and fields are
  #  0. start in ref scaffold
  #  1. end in ref scaffold
  #  2. start in query scaffold
  #  3. end in query scaffold
  #  4. match length in ref
  #  5. match length in query
  #  6. % identity of match
  #  7. total length of scaffold in ref
  #  8. total length of scaffold in query
  #  9. % coverage of contig in ref
  # 10. % coverage of contig in query
  # 11. scaffold UID in ref
  # 12. scaffold UID in query
    
}
close($ifh);









sub ReadScaffFile($)
{
  my $fname = shift;

  my %data;
  my $scaffUID = 0;
  my @fields;
  my $offset = 0;
  my $bps = 0;
  my $fauxKey;
  
  my $fh = new FileHandle $fname, "r"
    or die "Failed to open $fname for reading";
  while(<$fh>)
  {
    s/[\n\r\cZ]//g;

    @fields = split " ";
    
    # scaffold line starts with ">"
    if(substr($_,0,1) eq ">")
    {
      if($scaffUID != 0)
      {
        $fauxKey = "l" . "$scaffUID";
        $data{$fauxKey} = $offset - 100;
        $fauxKey = "b" . "$scaffUID";
        $data{$fauxKey} = $bps;
      }
      $scaffUID = substr($fields[0],1);
      $offset = 0;
      $bps = 0;
      next;
    }

    # contig line
    my @keepFields;
    # left end
    $keepFields[0] = $offset;
    # right end
    $keepFields[1] = $fields[2];
    # orientation
    $keepFields[2] = (($fields[1] eq "BE") ? 0 : 1);

    push @{$data{$scaffUID}}, [@keepFields];

    $offset += $fields[2] + $INTRA_SCAFFOLD_GAP_LENGTH;
    $bps += $fields[2];
  }
  close($fh);

  # store the last scaffolds length & numBPs
  $fauxKey = "l" . "$scaffUID";
  $data{$fauxKey} = $offset - 100;
  $fauxKey = "b" . "$scaffUID";
  $data{$fauxKey} = $bps;
  
  return %data;
}


sub DumpScaffData(%)
{
  my %data = @_;
  
  foreach my $scaffUID (sort {$a <=> $b} (keys(%data)))
  {
    next if(substr($scaffUID, 0, 1) eq "l" ||
            substr($scaffUID, 0, 1) eq "b");
    
    printf(">%s %d %d %d\n",
           $scaffUID,
           1 + $#{@{$data{$scaffUID}}},
           $data{"b" . "$scaffUID"},
           $data{"l" . "$scaffUID"});
    my $i;
    for($i = 0; $i <= $#{@{$data{$scaffUID}}}; $i++)
    {
      printf("%d %d %s\n",
             $data{$scaffUID}[$i][0],
             $data{$scaffUID}[$i][1],
             ($data{$scaffUID}[$i][2] == 0) ? "BE" : "EB");
    }
  }
}
