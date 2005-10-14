#!/usr/local/bin/perl
# $Id: runTampa.pl,v 1.4 2005-10-14 17:48:36 catmandew Exp $
#
# Wrapper to run and post-process results from TAMPA
# (Tool for Analyzing Mate Pairs in Assemblies)
#
# Ian M. Dew, Brian Walenz, Granger Sutton. "A Tool for Analyzing Mate
#   Pairs in Assemblies (TAMPA)", J Comp Biol 2005 12(5):497-513.
#
# Written by Ian Dew
#

use Carp;
# use strict;
use FileHandle;
use Getopt::Long;

# mate pair types & statuses
my %IntraIndices = ("mate" => 0,
                    "coincident" => 1,
                    "raw" => {"satisfied" => 2,
                              "stretched" => 3,
                              "compressed" => 4,
                              "outtie" => 5,
                              "normal" => 6,
                              "antinormal" => 7,},
                    "confirmed" => {"stretched" => 8,
                                    "compressed" => 9,
                                    "outtie" => -1,
                                    "normal" => -1,
                                    "antinormal" => -1,
                                    "inversion" => 10,
                                    "transposition" => 11,
                                    "insertions" => 23,
                                    "deletions" => 24,},
                    "below-threshold" => {"stretched" => 12,
                                          "compressed" => 13,
                                          "outtie" => -1,
                                          "normal" => -1,
                                          "antinormal" => -1,
                                          "inversion" => 14,
                                          "transposition" => 15,},
                    "polymorphic" => {"insertion" => 16,
                                      "inversion" => 17,},
                    "problematic" => {"stretched" => 18,
                                      "compressed" => 19,
                                      "outtie" => 20,
                                      "normal" => 21,
                                      "antinormal" => 22,
                                      "inversion" => -1,},);
my $LastIntraIndex = 24;

my %InterIndices = ("mate" => 0,
                    "stretched" => 1,
                    "outtie" => 2,
                    "normal" => 3,
                    "antinormal" => 4,);
my $LastInterIndex = 4;

my $lengthAccumulator = 0;
my %libMateCounter;

my %PARAMETERS = ("assemblyPrefix" => "",
                  "libFilename" => "",
                  "numSigmas" => 3,
                  "minPairs" => 2,
                  "binariesPath" => "",
                  "dontReestimate" => 0,
                  "reestIters" => 4,
                  "reestSigmas" => 4,
                  "dontDoIntra" => 0,
                  "dontDoInter" => 0,
                  "verboseLevel" => 0,
                  "gnuplotOutput" => 0,
                  "ataOutput" => 0,
                  "rawOutput" => 0,);

my $MY_VERSION = " Version 1.01 (Build " . (qw/$Revision: 1.4 $/ )[1]. ")";
my $MY_APPLICATION = "TAMPA";

my $REFERENCE = qq~
  Ian M. Dew, Brian Walenz, Granger Sutton, A Tool for Analyzing Mate
    Pairs in Assemblies (TAMPA), J Comp Bio. 2005 June; 12 (5):497-513.
  
~;

my $HELPTEXT = qq~
Run TAMPA (Tool for Analyzing Mate Pairs in Assemblies) on a genomic assembly.

    runTampa  [options]  <-a assembly>  <-l library>

    <-a assembly>  The prefix of the input filenames. Intra-sequence files
                   must have the form
                            assembly_([0-9]*)_intra.txt
                   Inter-sequence files must have the form
                            assembly_([0-9]*)_inter.txt
                   Refer to the user's manual for file formats.

    <-l library>  The file listing clone libraries, means, and standard
                  deviations.
                  Please refer to the user's manual for file formats.
  
    options:
      -h               Print help.
  
      -c               Print TAMPA citation.

      -v <level>       Set verbosity to level.

      -b <path>        Path to binaries.
                       Default = the current path.

      -r               Do NOT reestimate clone library means stddevs.
                       Conflicts with -i and -e.

      -i <iterations>  Iterate iterations times when reestimating clone
                       library means and stddevs.
                       Conflicts with -r.
                       DEFAULT = $PARAMETERS{"reestIters"}.

      -e <num_sigmas>  Exclude mate pairs beyond num_sigmas of the working
                       mean when reestimating library means and stddevs.
                       Conflicts with -r.
                       DEFAULT = $PARAMETERS{"reestSigmas"}.

      -s <num_sigmas>  To be satisfied, a mate pair must be within num_sigmas
                       of the library mean.
                       DEFAULT = $PARAMETERS{"numSigmas"}.

      -p <pairs>       The minimum number of mate pairs that must agree
                       to confirm an assembly problem.
                       DEFAULT = $PARAMETERS{"minPairs"}.

      -o               Do not run TAMPA on intra-sequence mate pairs.

      -x               Do not run TAMPA on inter-sequence mate pairs.

      -g               Generate output for gnuplot.

      -t               Generate ATAC output.

      -m               Generate raw output.

$MY_VERSION

~;


######################################################################
# Parse the command line
######################################################################
my $helpRequested;
my $printCitation;
my $reestIters = 0;
my $reestSigmas = 0;

GetOptions("a=s" => \$PARAMETERS{"assemblyPrefix"},
           "l=s" => \$PARAMETERS{"libFilename"},
           "h|help" => \$helpRequested,
           "c" => \$printCitation,
           "b=s" => \$PARAMETERS{"binariesPath"},
           "v|V|verbose:1" => \$PARAMETERS{"verboseLevel"},
           "r+" => \$PARAMETERS{"dontReestimate"},
           "i=i" => \$reestIters,
           "e=f" => \$reestSigmas,
           "s=f" => \$PARAMETERS{"numSigmas"},
           "p=i" => \$PARAMETERS{"minPairs"},
           "o+" => \$PARAMETERS{"dontDoIntra"},
           "x+" => \$PARAMETERS{"dontDoInter"},
           "g+" => \$PARAMETERS{"gnuplotOutput"},
           "t+" => \$PARAMETERS{"ataOutput"},
           "m+" => \$PARAMETERS{"rawOutput"}) or die $HELPTEXT;

if($helpRequested)
{
  print STDERR "Help requested:\n\n";
  print STDERR $HELPTEXT;
  exit 0;
}

if($printCitation)
{
  print STDOUT $REFERENCE;
  exit 0;
}

if(!$PARAMETERS{"assemblyPrefix"})
{
  print STDERR "Please specify an assembly\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if(!$PARAMETERS{"libFilename"})
{
  print STDERR "Please specify a library filename\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if($PARAMETERS{"dontReestimate"} && ($reestIters != 0 || $reestSigmas != 0))
{
  print STDERR "Reestimation spec(s) conflict with directive to not reestimate.\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if($PARAMETERS{"numSigmas"} <= 1)
{
  print STDERR "The number of sigmas must be positive\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if($PARAMETERS{"minPairs"} < 1)
{
  print STDERR "The minimum number of pairs to confirm must be positive\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if($PARAMETERS{"dontDoIntra"} && $PARAMETERS{"dontDoInter"})
{
  print STDERR "Processing neither intra- nor inter-sequence mate pairs means doing nothing.\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

if($PARAMETERS{"verboseLevel"})
{
  print "Assembly prefix:      " . $PARAMETERS{"assemblyPrefix"} . "\n";
  print "Library filename:     " . $PARAMETERS{"libFilename"} . "\n";
  print "Satisfieds are < +/-  " . $PARAMETERS{"numSigmas"} . " stddevs\n";
  print "Min pairs to confirm: " . $PARAMETERS{"minPairs"}. "\n";
  print "\n";

  if($PARAMETERS{"dontReestimate"})
  {
    print "Not reestimating clone library means and stddevs.\n";
  }
  else
  {
    print "Reestimating clone libraries in " . $PARAMETERS{"reestIters"} .
      " iterations using mate pairs within " . $PARAMETERS{"reestSigmas"} . " sigmas.\n";
  }
  print "\n";

  print "Running TAMPA on intra-sequence mate pairs.\n" if($PARAMETERS{"dontDoIntra"} == 0);
  print "Running TAMPA on inter-sequence mate pairs.\n" if($PARAMETERS{"dontDoInter"} == 0);
  print "\n";

  print "Binaries located in ";
  if($PARAMETERS{"binariesPath"})
  {
    print $PARAMETERS{"binariesPath"} . "\n";
  }
  else
  {
    print "current path.\n";
  }
  print "\n";
  
  print "Started ";
  system("date");
}

$PARAMETERS{"binariesPath"} .= "/" if($PARAMETERS{"binariesPath"});

######################################################################
# Get the complete list of sequences to process
######################################################################
my $currDir = ".";
opendir(DIR, $currDir) || die "can't opendir $currDir: $!";
my @entries = grep(-dr, readdir(DIR));
closedir DIR;

my %seqs;
my %intras;
my %inters;
foreach my $entry (@entries)
{
  # parse for intra and inter filenames
  next if($entry !~ /^($PARAMETERS{"assemblyPrefix"})_([a-zA-Z0-9]+)_int((ra)|(er))\.txt$/);
  my $seqNum = $2;
  $seqs{$seqNum} = 1;
  print "Found $entry\n" if($PARAMETERS{"verboseLevel"} > 1);
  
  $intras{$seqNum} = 1 if($entry =~ /intra/);
  $inters{$seqNum} = 1 if($entry =~ /inter/);
}

######################################################################
# Reestimate clone library means & stddevs
######################################################################
if(!$PARAMETERS{"dontReestimate"})
{
  my $rlfn = $PARAMETERS{"libFilename"} . ".Reest";
  print STDERR "Reestimating libraries and creating $rlfn.\n";

  # iterate through lines in library file
  my $numLibs = 0;
  my $ilfh = new FileHandle $PARAMETERS{"libFilename"}, "r" or
    die "Failed to open " . $PARAMETERS{"libFilename"} . " for reading";
  while(<$ilfh>)
  {
    s/[\n\r\cZ]//g;

    # Fields: ID, mean, stddev
    my @libFields = split;

    print "Reestimating library $libFields[0]\n"
      if($PARAMETERS{"verboseLevel"});
    
    # create a lengths file, populated with innie clone lengths
    my $ofn = $libFields[0] . "Lengths.txt";
    my $ofh = new FileHandle $ofn, "w" or
      die "Failed to open $ofn for writing";
    
    # iterate through all intra-sequence files to get innie clone lengths
    foreach my $intra (sort {$a <=> $b} (keys(%intras)))
    {
      my $iifn = $PARAMETERS{"assemblyPrefix"} . "_" . $intra . "_intra.txt";
      my $iifh = new FileHandle $iifn, "r" or
        die "Failed to open $iifn for reading";
      while(<$iifh>)
      {
        my @cloneFields = split;
        if($cloneFields[0] eq "I" && $cloneFields[3] == $libFields[0])
        {
          my $length = abs($cloneFields[5] - $cloneFields[4]);
          print $ofh "$length\n";
        }
      }
      close($iifh);
    }
    close($ofh);

    my $appender = ($numLibs == 0 ? " > " : " >> ");
    my $command = $PARAMETERS{"binariesPath"} . "reestimateLibs" .
      " -f " . $ofn .
      " -l " . $PARAMETERS{"libFilename"} .
      " -n " . $libFields[0] .
      " -s " . $PARAMETERS{"reestSigmas"} .
      " -i " . $PARAMETERS{"reestIters"} .
      $appender . $rlfn;
    
    print "Running $command\n" if($PARAMETERS{"verboseLevel"} > 1);
    system($command) == 0
      or die "Failed to run command\n$command\n\n";
    $numLibs++;
  }
  close($ilfh);

  $PARAMETERS{"libFilename"} = $rlfn;
}


######################################################################
# Run processIntra on each intra file
######################################################################
if(!$PARAMETERS{"dontDoIntra"})
{
  print "\n\nProcessing intra-sequence mate pairs\n\n"
    if($PARAMETERS{"verboseLevel"} > 1);

  # create concatenated breakpoints file
  my $bfn = $PARAMETERS{"assemblyPrefix"} . ".intra.breakpoints.tampa";
  my $bfh = new FileHandle $bfn, "w" or
      die "Failed to open $bfn for writing";
  
  # create csv file for sequence-by-sequence summary
  my $sfn = $PARAMETERS{"assemblyPrefix"} . ".intra.summary.tampa";
  my $sfh = new FileHandle $sfn, "w" or
      die "Failed to open $sfn for writing";

  # write title
  print $sfh "TAMPA Intra-sequence results for " .
    $PARAMETERS{"assemblyPrefix"} . " assembly\n\n";
  
  # write field labels
  print $sfh "\t" .
    "Excluded Clones\t\t" .
    "Raw Clones\t\t\t\t\t\t" .
    "Confirmed Problems\t\t\t\t" .
    "Below Threshold\t\t\t\t" .
    "Polymorphic\t\t" .
    "Problematic\t\t\t\t\t" .
    "Double Counted in Transpositions\t\n";

  # 25 fields after sequence id
  print $sfh "Sequence ID\t" .
    #    0            1
    "From bad lib\tCoincident\t" .
    #    2          3          4         5       6        7
    "Satisfied\tStretched\tCompressed\tOuttie\tNormal\tAntinormal\t" .
    #    8           9         10            11
    "Insertions\tDeletions\tInversions\tTranspositions\t" .
    #   12          13         14            15
    "Insertions\tDeletions\tInversions\tTranspositions\t" .
    #   16          17
    "Insertions\tInversions\t" .
    #   18         19        20      21        22
    "Stretched\tCompressed\tOuttie\tNormal\tAntinormal\t" .
    #   23          24
    "Insertions\tDeletions\n";

  # set up accumulators
  my @totals;
  my $i;
  for($i = 0; $i <= $LastIntraIndex; $i++)
  {
    $totals[$i] = 0;
  }
  
  foreach my $intra (sort {$a <=> $b} (keys(%intras)))
  {
    # run TAMPA on intra-sequence
    my $ofn = $PARAMETERS{"assemblyPrefix"} . "." . $intra . ".intra.summary.txt";
    my $ifn = $PARAMETERS{"assemblyPrefix"} . "_" . $intra . "_intra.txt";
    my $command = $PARAMETERS{"binariesPath"} . "processIntra" .
      " -l " . $PARAMETERS{"libFilename"} .
      " -m $ifn" .
      " -n " . $PARAMETERS{"numSigmas"} .
      " -f " . $PARAMETERS{"minPairs"};
    $command .= " -g " if($PARAMETERS{"gnuplotOutput"});
    $command .= " -a " if($PARAMETERS{"ataOutput"});
    $command .= " -r " if($PARAMETERS{"rawOutput"});
    $command .= " > " . $ofn;
    print "Running $command\n" if($PARAMETERS{"verboseLevel"} > 1);
    system($command) == 0
      or die "\nFailed to run command\n$command\n\n";

    # add breakpoints to concatenated breakpoints file
    $ifn = $PARAMETERS{"assemblyPrefix"} . "." . $intra . ".intra.breakpoints.txt";
    my $ifh = new FileHandle $ifn, "r" or
      die "Failed to open $ofn for reading";
    while(<$ifh>)
    {
      s/[\n\r\cZ]//g;
      printf($bfh "%s\t%s\n", $intra, $_);
    }
    
    # parse summary file & add to summary spreadsheet file
    printf $sfh "$intra";
    my @vals;
    for($i = 0; $i <= $LastIntraIndex; $i++)
    {
      $vals[$i] = 0;
    }
    
    my $iofh = new FileHandle $ofn, "r" or
      die "Failed to open $ofn for reading";
    while(<$iofh>)
    {
      my @fields = split;
      my $index = -1;
      next if($fields[1] eq "unsatisfied");

      if($fields[1] eq "raw" ||
         $fields[1] eq "confirmed" ||
         $fields[1] eq "below-threshold" ||
         $fields[1] eq "polymorphic" ||
         $fields[1] eq "problematic")
      {
        $index = $IntraIndices{$fields[1]}{$fields[2]};
      }
      elsif($fields[1] eq "is")
      {
        $lengthAccumulator += $fields[0];
        next;
      }
      elsif($fields[1] eq "clones")
      {
        $libMateCounter{$fields[4]} += $fields[0];
        next;
      }
      elsif(defined($IntraIndices{$fields[1]}))
      {
        $index = $IntraIndices{$fields[1]};
      }

      next if($index <= -1);
      $vals[$index] = $fields[0];
      $totals[$index] += $fields[0];
    }
    close($iofh);

    for($i = 0; $i <= $LastIntraIndex; $i++)
    {
      printf $sfh "\t$vals[$i]";
    }
    printf $sfh "\n";
  }
  close($bfh);

  # print totals
  printf $sfh "\nTotals";
  for($i = 0; $i <= $LastIntraIndex; $i++)
  {
    printf $sfh "\t$totals[$i]";
  }

  # calculate probability of detections & print
  
  printf $sfh "\n\nPlease cite the following in any publications" .
    "$REFERENCE\n";
  
  close($sfh);
}

if(!$PARAMETERS{"dontDoInter"})
{
  print "\n\nProcessing inter-sequence mate pairs\n\n"
    if($PARAMETERS{"verboseLevel"} > 1);

  # create concatenated breakpoints file
  my $bfn = $PARAMETERS{"assemblyPrefix"} . ".inter.breakpoints.tampa";
  my $bfh = new FileHandle $bfn, "w" or
      die "Failed to open $bfn for writing";
  
  # create csv file for sequence-by-sequence summary
  my $sfn = $PARAMETERS{"assemblyPrefix"} . ".inter.summary.tampa";
  my $sfh = new FileHandle $sfn, "w" or
      die "Failed to open $sfn for writing";

  # write title
  print $sfh "TAMPA Inter-sequence results for " .
    $PARAMETERS{"assemblyPrefix"} . " assembly\n\n";
  
  # write field labels
  print $sfh "\t\t" .
    "Confirmed intervals\n";

  # 5 fields after sequence id
  print $sfh "Sequence ID\t" .
    #    0
    "Raw clones\t" .
    #  1      2       3         4
    "Innie\tOuttie\tNormal\tAntinormal\n";

  # set up accumulators
  my @totals;
  my $i;
  for($i = 0; $i <= $LastInterIndex; $i++)
  {
    $totals[$i] = 0;
  }

  foreach my $inter (sort {$a <=> $b} (keys(%inters)))
  {
    my $ofn = $PARAMETERS{"assemblyPrefix"} . "." . $inter . ".inter.summary.txt";
    my $ifn = $PARAMETERS{"assemblyPrefix"} . "_" . $inter . "_inter.txt";
    my $command = $PARAMETERS{"binariesPath"} . "processInterCG" .
      " -l " . $PARAMETERS{"libFilename"} .
      " -m $ifn" .
      " -n " . $PARAMETERS{"numSigmas"} .
      " -f " . $PARAMETERS{"minPairs"};
    $command .= " -g " if($PARAMETERS{"gnuplotOutput"});
    $command .= " -a " if($PARAMETERS{"ataOutput"});
    $command .= " > " . $ofn;
    print "Running $command\n" if($PARAMETERS{"verboseLevel"} > 1);
    system($command) == 0
      or die "\nFailed to run command\n$command\n\n";

    # add breakpoints to concatenated breakpoints file
    $ifn = $PARAMETERS{"assemblyPrefix"} . "." . $inter . ".inter.breakpoints.txt";
    my $ifh = new FileHandle $ifn, "r" or
      die "Failed to open $ofn for reading";
    while(<$ifh>)
    {
      s/[\n\r\cZ]//g;
      printf($bfh "%s\t%s\n", $inter, $_);
    }
    
    # parse summary file & add to summary spreadsheet file
    printf $sfh "$inter";
    my @vals;
    for($i = 0; $i <= $LastInterIndex; $i++)
    {
      $vals[$i] = 0;
    }
    
    my $iofh = new FileHandle $ofn, "r" or
      die "Failed to open $ofn for reading";
    while(<$iofh>)
    {
      my @fields = split;

      $vals[$InterIndices{$fields[2]}] = $fields[0];
      $totals[$InterIndices{$fields[2]}] += $fields[0];
    }
    close($iofh);

    for($i = 0; $i <= $LastInterIndex; $i++)
    {
      printf $sfh "\t$vals[$i]";
    }
    printf $sfh "\n";
  }
  close($bfh);

  # print totals
  printf $sfh "\nTotals";
  for($i = 0; $i <= $LastInterIndex; $i++)
  {
    printf $sfh "\t$totals[$i]";
  }

  # calculate probability of detections & print
  
  printf $sfh "\n\nPlease cite the following in any publications" .
    "$REFERENCE\n";
  
  close($sfh);
}


if($PARAMETERS{"verboseLevel"})
{
  print "Done ";
  system("date");
}
