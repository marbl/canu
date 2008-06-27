#!/usr/bin/env perl
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
# $Id: relabelCloneLibs.pl,v 1.4 2008-06-27 06:29:17 brianwalenz Exp $
#

use FileHandle;
use Getopt::Long;

# user specifies
#  lib of clones (left UID, right UID) to relabel
#  lib to relabel them as
#  file in which they are mislabelled

# read file of clones lengths that should be lib a instead of lib b
#   into hashtable
# read file with mislabelled clones
#   change library as appropriate
#   write to stdout

sub usage
{
  printf STDERR "Usage: ${0}  [-h]  -c clonesFile  -l libNum  files...\n" .
    "  -h               print usage\n" .
    "  -c clonesFile    name of file listing mislabelled clones\n" .
    "                     format is length  lowFragUID  highFragUID\n" .
    "  -l libNum        UID of library the clones should belong to\n" .
    "  files...         list of TAMPA files in which the clones are mislabelled\n\n";
  exit;
}

my $clonesFilename = "";
my $libNum = 0;
my $help = 0;

GetOptions("c=s", => \$clonesFilename,
           "l=i", => \$libNum,
           "h", => \$help) or die "Option error";

usage if($help != 0 || $clonesFilename eq "" || $libNum == 0);

# hashtable:
# key = lowFragID
# value = whatever
my %mislabelled;

# read in list of mislabelled clones
my $fh;
$fh = new FileHandle $clonesFilename, "r" or die "Failed to open $clonesFilename for reading";
while(<$fh>)
{
  s/[\n\r\cZ]//g;

  my @fields = split " ";
  $mislabelled{$fields[1]} = $fields[2];
}
close($fh);

# read tampa files and modify them
my $i;
for($i = 0; $i <= $#ARGV; $i++)
{
  my $ifn = $ARGV[$i];
  my $ifh = new FileHandle $ifn, "r" or die "Failed to open $ifn for reading";

  my $ofn = "/tmp/$ARGV[$i]";
  my $ofh = new FileHandle $ofn, "w" or die "Failed to open $ofn for writing";

  my $numChanged = 0;

  # different format for intra vs inter files
  my $type = (($ARGV[$i] =~ m/_intra.txt$/) ? 0 : 1);

  while(<$ifh>)
  {
    s/[\n\r\cZ]//g;

    my @fields = split " ";
    if($type == 0)
    {
      if(defined($mislabelled{$fields[1]}) ||
         defined($mislabelled{$fields[2]}))
      {
          $fields[3] = $libNum;
          $numChanged++;
      }
    }
    else
    {
      if(defined($mislabelled{$fields[0]}) ||
         defined($mislabelled{$fields[4]}))
      {
          $fields[8] = $libNum;
          $numChanged++;
      }
    }

    printf($ofh "%s", $fields[0]);
    my $j;
    for($j = 1; $j <= $#fields; $j++)
    {
      printf($ofh " %s", $fields[$j]);
    }
    printf($ofh "\n");
  }
  close($ifh);
  close($ofh);

  printf(STDERR "Changed $numChanged in $ifn\n");

  my $command = "mv -f $ofn $ifn";
  if(system($command))
  {
    printf(STDERR "Failed to run $command\n");
    exit;
  }
}
