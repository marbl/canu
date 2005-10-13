#!/usr/local/bin/perl
# $Id: assemblyCompare.pl,v 1.2 2005-10-13 21:21:12 catmandew Exp $
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

my $MY_VERSION = " Version 1.01 (Build " . (qw/$Revision: 1.2 $/ )[1]. ")";
my $MY_APPLICATION = "assemblyCompare";

my $HELPTEXT = qq~
Compare two assemblies

    assemblyCompare  [options]  -d refDir  -a refAssembly
                                -d queryDir  -a queryAssembly

    -d refDir        The 'reference' assembly directory

    -a refAssembly   The 'reference' assembly name

    -d queryDir        The 'query' assembly directory

    -a queryAssembly   The 'query' assembly name

    Reference directory and assembly name must both be listed before
    the query reference directory and assembly name.
  
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

if($#dirs != 1 || ! -d $dirs[0] || ! -d $dirs[1])
{
  print STDERR "Please specify a reference and query directory\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if($#assemblies != 1)
{
  print STDERR "Please specify a reference and query assembly name\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

printf("\nComparing assembly %s in %s\n", $assemblies[0], $dirs[0]);
printf("with assembly %s in %s\n\n", $assemblies[1], $dirs[1]);


######################################################################
# Compare qc files
######################################################################
if(!$dontQC)
{
  my $rfn = $dirs[0] . "/" . $assemblies[0] . ".qc";
  my $qfn = $dirs[1] . "/" . $assemblies[1] . ".qc";

  if(-f $rfn && -f $qfn)
  {
    printf("==========> QC Statistics Comparison\n\n");
    my $command = "qcCompare.pl -f $rfn -f $qfn";
    printf STDERR "Running $command\n";
    system($command);
  }    
}


if(!$dontTampa)
{
  my @types = ("intra", "inter");
  # assume either results files are present or TAMPA has not been run
  my $tampRun = 1;
  for($i = 0; $i <= $#dirs; $i++)
  {
    for($j = 0; $j <= $#types; $j++)
    {
      my $fn = $dirs[$i] . "/" . $assemblies[$i] . "." .
        $types[$j] . ".summary.tampa";
      if(! -f $fn )
      {
        $tampaRun = 0;
        last;
      }
      last if($tampaRun == 0);
    }
  }

  if(!$tampaRun)
  {
    # run tampa
    for($i = 0; $i <= $#dirs; $i++)
    {
      chdir($dirs[$i]);
      my $command = "asm2TampaResults -a $assemblies[$i]";
      printf STDERR "Running $command\n";
      system($command);
    }
    chdir($PWD);
  }

  my $command = "tampaCompare.pl -d $dirs[0] -a $assemblies[0] -d $dirs[1] -a $assemblies[1]";
  printf STDERR "Running $command\n";
  system($command);
}


if(!dontMummer)
{
  # assume mummer hasn't been run on this pair, since a pairwise thing..
  my $pid = getppid();
  my $prefix = "$dirs[0]_$dirs[1]_$pid";
  my @fnps;
  for(my $i = 0; $i <= $#dirs; $i++)
  {
    $fnps[$i] = $dirs[$i] . "/" . $assemblies[i]"
  }
  my $command = "nucmer -p $prefix $fnps[0].scaffolds.fasa $fnps[1].scaffolds.fasta";
  printf STDERR "Running $command\n";
  system($command);

  # now generate secondary output
  $command = "show-coords -THcl -I 99 $prefix.delta > $prefix.show-coords";
  printf STDERR "Running $command\n";
  system($command);

  # now compare
  $command = "mummerCompare.pl -r $fnps[0].scaff -q $fnps[1].scaff -c $prefix.show-coords";
  printf STDERR "Running $command\n";
  system($command);
}
