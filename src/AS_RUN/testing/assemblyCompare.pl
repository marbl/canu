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
# $Id: assemblyCompare.pl,v 1.6 2005-12-16 22:12:38 catmandew Exp $
#

# Program to compare two assemblies
#
#  Compare two qc files
#  Compare TAMPA results
#  Run nucmer & show-coords & show-snps
#  Compare show-coords & show-snps output
#
#   Written by Ian Dew
#

use Carp;
use strict;
use FileHandle;
use Getopt::Long;
use Env qw(PWD);

my $MY_VERSION = " Version 1.01 (Build " . (qw/$Revision: 1.6 $/ )[1]. ")";
my $MY_APPLICATION = "assemblyCompare";

my $HELPTEXT = qq~
Compare two assemblies

    assemblyCompare  [options]  -d refDir  -a refAssembly
                                -d queryDir  -a queryAssembly
                                ...

    -d refDir        The 'reference' assembly directory

    -a refAssembly   The 'reference' assembly name

    -d queryDir        The 'query' assembly directory

    -a queryAssembly   The 'query' assembly name

    Reference directory and assembly name must both be listed before
    the query reference directories and assembly names. Any number of
    query directories/assemblies may be specified.
  
    options:
      -h               Print help.
  
      -v <level>       Set verbosity to level.

      -c               Do not compare qc files

      -t               Do not compare TAMPA results

      -m               Do not run mummer and compare results

$MY_VERSION

~;


######################################################################
# Parse the command line
######################################################################
my $helpRequested;
my $verboseLevel = 0;
my @dirs;
my @assemblies;
my $dontQC;
my $dontTampa;
my $dontMummer;

GetOptions("d=s" => \@dirs,
           "a=s" => \@assemblies,
           "c" => \$dontQC,
           "t" => \$dontTampa,
           "m" => \$dontMummer,
           "h|help" => \$helpRequested,
           "v|V|verbose:1" => \$verboseLevel
           ) or die $HELPTEXT;

if($helpRequested)
{
  print STDERR "Help requested:\n\n";
  print STDERR $HELPTEXT;
  exit 0;
}

if($#dirs < 1 || ! -d $dirs[0] || ! -d $dirs[1])
{
  print STDERR "Please specify a reference and query directory\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if($#assemblies < 1 || $#assemblies != $#dirs)
{
  print STDERR "Please specify at least one reference and query assembly name,\n";
  print STDERR "and please specify the same number of dirs as assemblies\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}


######################################################################
# Compare qc files
######################################################################
if(!$dontQC)
{
  my $rfn = $dirs[0] . "/" . $assemblies[0] . ".qc";
  for(my $i = 1; $i <= $#dirs; $i++)
  {
    my $qfn = $dirs[$i] . "/" . $assemblies[$i] . ".qc";

    if(-f $rfn && -f $qfn)
    {
      printf("==========> QC Statistics Comparison: Reference vs Assembly $i\n\n");
      my $command = "qcCompare -f $rfn -f $qfn";
      printf STDERR "Running $command\n";
      my $retval = system($command);
      if($retval != 0)
      {
        printf STDERR "Command failed. Aborting.\n";
        exit 1;
      }
    }
  }
}


if(!$dontTampa)
{
  my @types = ("intra", "inter");
  # assume either results files are present or TAMPA has not been run
  for(my $i = 0; $i <= $#dirs; $i++)
  {
    my $generateTAMPAData = 0;
    for(my $j = 0; $j <= $#types; $j++)
    {
      my $fn = $dirs[$i] . "/" . $assemblies[$i] . "." .
        $types[$j] . ".summary.tampa";
      if(! -f $fn )
      {
        $generateTAMPAData = 1;
        last;
      }
      last if($generateTAMPAData == 1);
    }

    if($generateTAMPAData)
    {
      # run tampa
      chdir($dirs[$i]);
      my $command = "asm2TampaResults -a $assemblies[$i]";
      printf STDERR "Running $command\n";
      my $retval = system($command);
      if($retval != 0)
      {
        printf STDERR "Command failed. Aborting.\n";
        exit 1;
      }

    }
    chdir($PWD);
  }

  for(my $i = 1; $i <= $#dirs; $i++)
  {
    printf("==========> TAMPA Comparison: Reference vs. Assembly $i\n\n");
    my $command = "tampaCompare -d $dirs[0] -a $assemblies[0] -d $dirs[$i] -a $assemblies[$i]";
    printf STDERR "Running $command\n";
    my $retval = system($command);
    if($retval != 0)
    {
      printf STDERR "Command failed. Aborting.\n";
      exit 1;
    }
  }
}

#$dontMummer = 1;
if(!$dontMummer)
{
  my $command;
  my @fnps;
  $fnps[0] = $dirs[0] . "/" . $assemblies[0];
  for(my $i = 1; $i <= $#dirs; $i++)
  {
    $fnps[$i] = $dirs[$i] . "/" . $assemblies[$i];
    
    my $prefix = $assemblies[0] . "_" . $assemblies[$i] . "_" . $i . "_nucmer";

    my $clusterFN = $prefix . ".cluster";
    my $deltaFN = $prefix . ".delta";
    if( ! -f $clusterFN || ! -f $deltaFN )
    {
      $command = "nucmer -p $prefix $fnps[0].scaffolds.fasta $fnps[$i].scaffolds.fasta";
      printf STDERR "Running $command\n";
      my $retval = system($command);
      if($retval != 0)
      {
        printf STDERR "Command failed. Aborting.\n";
        exit 1;
      }
    }
    
    my $showCoordsFN = $prefix . ".show-coords";
    if( ! -f $showCoordsFN )
    {
      # now generate secondary output
      $command = "show-coords -THcl -I 99 $prefix.delta > $prefix.show-coords";
      printf STDERR "Running $command\n";
      my $retval = system($command);
      if($retval != 0)
      {
        printf STDERR "Command failed. Aborting.\n";
        exit 1;
      }
    }
    
    # now compare
    printf("==========> Mummer Comparison: Reference vs. Assembly $i\n\n");
    $command = "analyzeMummerMapping -r $fnps[0].scaff -q $fnps[$i].scaff -s $showCoordsFN";
    printf STDERR "Running $command\n";
    my $retval = system($command);
    if($retval != 0)
    {
      printf STDERR "Command failed. Aborting.\n";
      exit 1;
    }
  }
}
