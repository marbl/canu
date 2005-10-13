#!/usr/local/bin/perl
# $Id: tampaCompare.pl,v 1.2 2005-10-13 21:21:12 catmandew Exp $
#
# Program to compare tampa output files
#
# Just compare totals for now
#
# Intra-summary totals line is tab-delimited with the following fields:
#  0. "Totals"
#  1. Excluded clones from bad libs
#  2. Excluded coincident clones
#  3. Raw satisfied
#  4. Raw stretched
#  5. Raw compressed
#  6. Raw outtie
#  7. Raw normal
#  8. Raw antinormal
#  9. Confirmed insertions
# 10. Confirmed deletions
# 11. Confirmed inversions
# 12. Confirmed tranpsositions
# 13. Below-threshold insertions
# 14. Below-threshold deletions
# 15. Below-threshold inversions
# 16. Below-threshold transpositions
# 17. Polymorphic insertions
# 18. Polymorphic inversions
# 19. Problematic stretched
# 20. Problematic compressed
# 21. Problematic outtie
# 22. Problematic normal
# 23. Problematic antinormal
# 24. Double-counted insertions in transpositions
# 25. Double-counted deletions in transpositions
#
# Most interesting are fields: 2-12, 17-18
#
#
# Inter-summary totals line is tab-delimited with the following fields:
#  0. "Totals"
#  1. Raw inter-sequence clones
#  2. Confirmed innie
#  3. Confirmed outtie
#  4. Confirmed normal
#  5. Confirmed antinormal
#
# All are of interest
#
#   Written by Ian Dew
#

use Carp;
use strict;
use FileHandle;
use Getopt::Long;

my $MY_VERSION = " Version 1.01 (Build " . (qw/$Revision: 1.2 $/ )[1]. ")";
my $MY_APPLICATION = "tampaCompare";

my $HELPTEXT = qq~
Compare two sets of TAMPA output files

    tampaCompare  [options]  -d refDir  -a refAssembly
                             -d queryDir1  -a queryAssembly1
                             ...

    -d refDir        The directory with the 'reference' genome
  
    -a refAssembly   The 'reference' assembly name prefix
                        These two files must exist:
                           refDir/refAssembly.intra.summary.tampa
                           refDir/refAssembly.inter.summary.tampa

    -d queryDir1        The directory with the first 'query' genome
  
    -a queryAssembly1   The first 'query' assembly name prefix
                        These two files must exist:
                           queryDir1/queryAssembly1.intra.summary.tampa
                           queryDir1/queryAssembly1.inter.summary.tampa


    An arbitrary number of query directories and assembly names
      may be specified. Directories and assembly names must be paired.

    options:
      -h               Print help.
  
      -v <level>       Set verbosity to level.

$MY_VERSION

~;

my @LTYPES = ("intra", "inter");
my @UTYPES = ("Intra", "Inter");

######################################################################
# Parse the command line
######################################################################
my $helpRequested;
my $verboseLevel = 0;
my @dirs;
my @assemblies;

GetOptions("d=s" => \@dirs,
           "a=s" => \@assemblies,
           "h|help" => \$helpRequested,
           "v|V|verbose:1" => \$verboseLevel
           ) or die $HELPTEXT;

if($helpRequested)
{
  print STDERR "Help requested:\n\n";
  print STDERR $HELPTEXT;
  exit 0;
}

# at least 2 assemblies and 
if($#dirs < 1 || $#dirs != $#assemblies)
{
  print STDERR "Please specify two or more directories & assembly names\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

# check that all required files exist
for(my $i = 0; $i <= $#dirs; $i++)
{
  for(my $j = 0; $j <= $#LTYPES; $j++)
  {
    my $fname = $dirs[$i] . "/" . $assemblies[$i] . "." .
      $LTYPES[$j] . ".summary.tampa";
    if(! -f $fname)
    {
      print STDERR "%s is not a file!\n", $fname;
      print STDERR $HELPTEXT;
      exit 1;
    }
  }
}


######################################################################
# Define fields are interesting
######################################################################
my @interest = ({"Excluded coincident clones" =>  2,
                 "Raw satisfied" =>  3,
                 "Raw stretched" =>  4,
                 "Raw compressed" =>  5,
                 "Raw outtie" =>  6,
                 "Raw normal" =>  7,
                 "Raw antinormal" =>  8,
                 "Confirmed insertions" =>  9,
                 "Confirmed deletions" => 10,
                 "Confirmed inversions" => 11,
                 "Confirmed tranpsositions" => 12,
                 "Polymorphic insertions" => 17,
                 "Polymorphic inversions" => 18,},
                {"Raw inter-sequence clones" =>  1,
                 "Confirmed innie" =>  2,
                 "Confirmed outtie" =>  3,
                 "Confirmed normal" =>  4,
                 "Confirmed antinormal" =>  5,});


######################################################################
# Read in the reference intra & inter files
######################################################################
my @rFields;
for(my $j = 0; $j <= $#LTYPES; $j++)
{
  @{$rFields[$j]} = GetTAMPATotals($dirs[0], $assemblies[0], $LTYPES[$j]);
}


######################################################################
# Loop over query directories & compare
######################################################################
my $header = "Category\tRef#\tQuery#\tDelta#\tDelta\%";
for(my $i = 1; $i <= $#dirs; $i++)
{
  printf "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
  printf("Reference Directory and Assembly: %s\t%s\n",
         $dirs[0], $assemblies[$0]);
  printf("Query Directory and Assembly:     %s\t%s\n",
         $dirs[$i], $assemblies[$i]);
  
  # loop over intra & inter
  for(my $j = 0; $j <= $#LTYPES; $j++)
  {
    printf("\n%s\n", $UTYPES[$j]);

    my $printedHeader = 0;
    my @qFields = GetTAMPATotals($dirs[$i], $assemblies[$i], $LTYPES[$j]);

    # loop over fields of interest
    foreach my $foi (sort(keys(%{$interest[$j]})))
    {
      if($rFields[$j][$interest[$j]{$foi}] !=
         $qFields[$interest[$j]{$foi}])
      {
        if($printedHeader == 0)
        {
          printf("%s\n", $header);
          $printedHeader = 1;
        }
        my $delta = $qFields[$interest[$j]{$foi}] -
          $rFields[$j][$interest[$j]{$foi}];
        my $pctDelta = 100;
        if($rFields[$j][$interest[$j]{$foi}] != 0)
        {
          $pctDelta = 100 * $delta / $rFields[$j][$interest[$j]{$foi}];
        }
        
        printf("%s\t%d\t%d\t%.2f\t%.2f\n",
               $foi,
               $rFields[$j][$interest[$j]{$foi}],
               $qFields[$interest[$j]{$foi}],
               $delta, $pctDelta);
      }
    }
  }
}


sub GetTAMPATotals($$$)
{
  my $dir = shift;
  my $name = shift;
  my $type = shift;
  my @fields;

  my $ifn = $dir . "/" . $name . "." . $type . ".summary.tampa";
  my $ifh = new FileHandle $ifn, "r"
    or die "Failed to open $ifn for reading";
  while(<$ifh>)
  {
    s/[\n\r\cZ]//g;
    
    @fields = split "\t";
    last if($fields[0] eq "Totals");
  }
  close($ifh);

  return @fields;
}
