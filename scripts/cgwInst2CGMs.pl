#!/usr/bin/perl -w
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

# Mates:
# Have intra-unitig, intra-contig, intra-scaffold, & totals
#    by celera reads, BAC ends, external reads, & totals
# 16 files with 4(x2) columns each

# unitig sizes, surrogate sizes

# scaffold size, gaps, contig sizes


# Script to convert instrumenter output with verbose3 option
# to a set of celagram files
# interested in actual numbers and ratios. output is by scaffold
# file0: all-mates
#  scfIID happy# miso# close# far# bad# inter-scaffold# mateless# happy% ....
# file1: intra-mates
#  scfIID happy# miso# close# far# bad# inter-contig# mateless# happy% ....
# file2: inter-mates
#  scfIID happy# miso# close# far# bad# inter-scaffold# mateless# happy% ....
# file3: contig data - sizes, unitigs, surrogates, fragments
#  scfIID #unitigs minU meanU maxU #surrogates minS meanS maxS
# file4: scaffold data - sizes, contigs, gaps
#  scfIID size #contigs minC meanC maxC #gaps #neg #pos minP meanP maxP
# file5: fragment data
#  scfIID #reads #BACends #extFrgs #locales

# file6: read mate data
#  scfIID #happy #miso #close #far #bad %happy ....
# file7: BAC end mate data
#  scfIID #happy #miso #close #far #bad %happy ....
# file8: external read mate data
#  scfIID #happy #miso #close #far #bad %happy ....
# file9: external fragment mate data
#  scfIID #happy #miso #close #far #bad %happy ....
# file10: mate type happiness data
#  scfIID fragHappy%  extHappy%  beHappy%  miso, close, far, bad

use Carp;
use strict 'vars';
use Getopt::Std;

use vars qw($opt_i $opt_o);

###########################################################
# forward declarations
sub CloseFiles;
sub WriteToAllFiles;
sub PrintToAllMatesFile;
sub PrintToInterFile;
sub PrintToScaffoldFile;
sub PrintToIntraFile;
sub PrintToContigFile;
sub PrintToFrgFile;
sub PrintToHappyMatesFile;
sub PrintToMisoMatesFile;
sub PrintToFarMatesFile;
sub PrintToCloseMatesFile;
sub PrintToBadMatesFile;
sub PrintToReadMatesFile;
sub PrintToExternalMatesFile;
sub PrintToBacEndMatesFile;
sub PrintToMateTypesFile;
sub ComputeMatePercentages;
###########################################################


###########################################################
# Check that command line was sufficient
my $Usage = "Usage: $0
-i   instrumenting text in
-o   output filename stub
";

getopts('i:o:');

if(! $opt_i || ! $opt_o)
{
  print $Usage;
  exit(0);
}
###########################################################


###########################################################
# Open input & output files
# define names
my $unitigMatesOut = "${opt_o}_unitigMates.cgm";
my $contigMatesOut = "${opt_o}_contigMates.cgm";
my $scaffoldMatesOut = "${opt_o}_scaffoldMates.cgm";
my $allMatesOut = "${opt_o}_allMates.cgm";
my $happyMatesOut = "${opt_o}_happyMates.cgm";
my $misoMatesOut = "${opt_o}_misoMates.cgm";
my $farMatesOut = "${opt_o}_farMates.cgm";
my $closeMatesOut = "${opt_o}_closeMates.cgm";
my $badMatesOut = "${opt_o}_badMates.cgm";
my $readMatesOut = "${opt_o}_readMates.cgm";
my $bacEndMatesOut = "${opt_o}_bacEndMates.cgm";
my $extReadMatesOut = "${opt_o}_extReadsMates.cgm";
my $extFragMatesOut = "${opt_o}_extFragsMates.cgm";
my $mateTypesOut = "${opt_o}_mateTypes.cgm";

my $ctgOut = "${opt_o}_contigs.cgm";
my $scfOut = "${opt_o}_scaffolds.cgm";
my $frgOut = "${opt_o}_fragments.cgm";

#open files
# input
open(TXTIN, "$opt_i") || die "Failed to open file $opt_i for reading\n";

# summary/totals
open(UMOUT, "> $unitigMatesOut") || die "Failed to open file $unitigMatesOut for writing\n";
open(CMOUT, "> $contigMatesOut") || die "Failed to open file $contigMatesOut for writing\n";
open(SMOUT, "> $scaffoldMatesOut") || die "Failed to open file $scaffoldMatesOut for writing\n";
open(AMOUT, "> $allMatesOut") || die "Failed to open file $allMatesOut for writing\n";

open(HAPPYOUT, "> $happyMatesOut") || die "Failed to open file $happyMatesOut for writing\n";
open(MISOOUT, "> $misoMatesOut") || die "Failed to open file $misoMatesOut for writing\n";
open(FAROUT, "> $farMatesOut") || die "Failed to open file $farMatesOut for writing\n";
open(CLOSEOUT, "> $closeMatesOut") || die "Failed to open file $closeMatesOut for writing\n";
open(BADOUT, "> $badMatesOut") || die "Failed to open file $badMatesOut for writing\n";

open(READOUT, "> $readMatesOut") || die "Failed to open file $readMatesOut for writing\n";
open(BEOUT, "> $bacEndMatesOut") || die "Failed to open file $bacEndMatesOut for writing\n";
open(EXTROUT, "> $extReadMatesOut") || die "Failed to open file $extReadMatesOut for writing\n";
open(EXTFOUT, "> $extFragMatesOut") || die "Failed to open file $extFragMatesOut for writing\n";
open(TYPEOUT, "> $mateTypesOut") || die "Failed to open file $mateTypesOut for writing\n";

# contig, scaffold stats, #s of fragments by type
open(CTGOUT, "> $ctgOut") || die "Failed to open file $ctgOut for writing\n";
open(SCFOUT, "> $scfOut") || die "Failed to open file $scfOut for writing\n";
open(FRGOUT, "> $frgOut") || die "Failed to open file $frgOut for writing\n";

# write celagram header lines
print UMOUT "scfIID happy# miso# close# far# bad# inter# mateless# happy%...\n";
print CMOUT "scfIID happy# miso# close# far# bad# inter# mateless# happy%...\n";
print SMOUT "scfIID happy# miso# close# far# bad# inter# mateless# happy%...\n";
print AMOUT "scfIID happy# miso# close# far# bad# inter# mateless# happy%...\n";

print HAPPYOUT "scfIID unitig# contig# scaffold# total# unitig%...\n";
print MISOOUT "scfIID unitig# contig# scaffold# total# unitig%...\n";
print FAROUT "scfIID unitig# contig# scaffold# total# unitig%...\n";
print CLOSEOUT "scfIID unitig# contig# scaffold# total# unitig%...\n";
print BADOUT "scfIID unitig# contig# scaffold# total# unitig%...\n";

print READOUT "scfIID happy# miso# close# far# happy%...\n";
print BEOUT "scfIID happy# miso# close# far# happy%...\n";
print EXTROUT "scfIID happy# miso# close# far# happy%...\n";
print EXTFOUT "scfIID happy# miso# close# far# happy%...\n";
print TYPEOUT "scfIID happy%read  happy%be  happy%extReads  happy%extFrags  miso%s  close%s  far%s  bad%s\n";

print CTGOUT "scfIID #unitigs minU meanU maxU #surrogates minS meanS maxS\n";
print SCFOUT "scfIID size #contigs minC meanC maxC #gaps #negG #posG minP meanP maxP\n";
print FRGOUT "scfIID #reads #BACends #externalReads #extFrgs #locales\n";
###########################################################


###########################################################
# variables
# for reading input lines
my $newline;

# for holding instrumenter values
my $scaffoldID = 0;

# fragment types with mates
my $allFrags = 4;

# unitig-mates
my @happyUnitig = (0) x ($allFrags + 1);
my @misoUnitig = (0) x ($allFrags + 1);
my @farUnitig = (0) x ($allFrags + 1);
my @closeUnitig = (0) x ($allFrags + 1);
my @matelessUnitig = (0) x ($allFrags + 1);
my $interUnitig = 0;

# contig-mates
my @happyContig = (0) x ($allFrags + 1);
my @misoContig = (0) x ($allFrags + 1);
my @farContig = (0) x ($allFrags + 1);
my @closeContig = (0) x ($allFrags + 1);
my @matelessContig = (0) x ($allFrags + 1);
my $interContig = 0;

# scaffold-mates
my @happyScaffold = (0) x ($allFrags + 1);
my @misoScaffold = (0) x ($allFrags + 1);
my @farScaffold = (0) x ($allFrags + 1);
my @closeScaffold = (0) x ($allFrags + 1);
my @matelessScaffold = (0) x ($allFrags + 1);
my $interScaffold = 0;

# contigs
my $numUnitigs = 0;
my $minUnitig = 0;
my $meanUnitig = 0;
my $maxUnitig = 0;
my $numSurrogates = 0;
my $minSurrogate = 0;
my $meanSurrogate = 0;
my $maxSurrogate = 0;
  
# scaffolds
my $size = 0;
my $numContigs = 0;
my $minContig = 0;
my $meanContig = 0;
my $maxContig = 0;
my $numGaps = 0;
my $numNegativeGaps = 0;
my $numPositiveGaps = 0;
my $minPositiveGap = 0;
my $meanPositiveGap = 0;
my $maxPositiveGap = 0;
  
# fragments
my $numReads = 0;
my $numBacEnds = 0;
my $numExtReads = 0;
my $numExtFrags = 0;
my $numLocales = 0;

# summary variables
my $numHappy = 0;
my $numMiso = 0;
my $numClose = 0;
my $numFar = 0;
my $numBad = 0;
my $numInter = 0;
my $numMateless = 0;
my $numMates = 0;
my $pctHappy = 0;
my $pctMiso = 0;
my $pctClose = 0;
my $pctFar = 0;
my $pctBad = 0;
my $pctInter = 0;

my $pctUnitig = 0;
my $pctContig = 0;
my $pctScaffold = 0;

# to track which block in parsing instrumenting file
my $inScaffold = 0;
my $inUnitig = 1;
my $inContig = 2;
my $whichBlock = $inScaffold;

my $inGapSizes = 0;
my $inUnitigSizes = 0;
my $inMateSummary = 0;

my $inReads = 0;
my $inBacEnds = 1;
my $inExtReads = 2;
my $inExtFrags = 3;
my $inTotal = 4;
my $inMateless = 5;
my $inExternal = 6;
my $whichFrags = $inReads;
###########################################################


###########################################################
# loop over instrumenting input file
while($newline = <TXTIN>)
{
  if($newline =~ /Unitig summary:/)
  {
    $whichBlock = $inUnitig;
    next;
  }
  elsif($newline =~ /Contig summary:/)
  {
    $whichBlock = $inContig;
    next;
  }
  elsif($newline =~ /Instrumenting Scaffold (\d+)/)
  {
    $scaffoldID = $1;
    $whichBlock = $inScaffold;
    next;
  }
  
  # if we are in a scaffold block (top and bottom)
  if($whichBlock == $inScaffold)
  {
#    print "In scaffold $scaffoldID\n";
    if($newline =~ /Size: (\d+)/)
    {
      $size = $1;
    }
    elsif($newline =~ /Scaffold gap sizes/)
    {
      $inGapSizes = 1;
    }
    elsif($newline =~ /Contig sizes/)
    {
      $inGapSizes = 0;
    }
    elsif($newline =~ /Celera reads:/)
    {
#      print "found Celera reads. inMateSummary = $inMateSummary\n";
      $inMateSummary = $inMateSummary + 1;
      $whichFrags = $inReads;
    }
    elsif($inMateSummary == 2)
    {
      # detect the block
#      print "In matesummary == 2\n";
      if($newline =~ /BAC ends:/)
      {
        $whichFrags = $inBacEnds;
        next;
      }
      elsif($newline =~ /External reads:/)
      {
        $whichFrags = $inExtReads;
        next;
      }
      elsif($newline =~ /External frags:/)
      {
        $whichFrags = $inExtFrags;
        next;
      }
      elsif($newline =~ /Total:/)
      {
        $whichFrags = $inTotal;
        next;
      }
      elsif($newline =~ /Mateless:/)
      {
        $whichFrags = $inMateless;
        next;
      }
      elsif($newline =~ /External mates:/)
      {
        $whichFrags = $inExternal;
        next;
      }

      # detect which line in the block
      if($whichFrags == $inMateless)
      {
        if($newline =~ /(\d+) reads/)
        {
          $matelessScaffold[$inReads] = $1;
        }
        elsif($newline =~ /(\d+) BAC ends/)
        {
          $matelessScaffold[$inBacEnds] = $1;
        }
        elsif($newline =~ /(\d+) external reads/)
        {
          $matelessScaffold[$inExtReads] = $1;
        }
        elsif($newline =~ /(\d+) external frags/)
        {
          $matelessScaffold[$inExtFrags] = $1;
        }
        elsif($newline =~ /(\d+) total/)
        {
          $matelessScaffold[$inTotal] = $1;
        }
      }
      elsif($whichFrags == $inExternal)
      {
        if($newline =~ /(\d+) fragments with external mates/)
        {
          $interScaffold = $1;
        }
        elsif($newline =~ /(\d+) locales/)
        {
          $numLocales = $1;
          
          # print all-mates, inter-mates, & scaffold data
          WriteToAllFiles();
          $inMateSummary = 0;
        }
      }
      else
      {
        if($newline =~ /(\d+) happy inter-contig/)
        {
          $happyScaffold[$whichFrags] = $1;
        }
        elsif($newline =~ /(\d+) mis-oriented inter-contig/)
        {
          $misoScaffold[$whichFrags] = $1;
        }
        elsif($newline =~ /(\d+) mis-separated too close inter-contig/)
        {
          $closeScaffold[$whichFrags] = $1;
        }
        elsif($newline =~ /(\d+) mis-separated too far inter-contig/)
        {
          $farScaffold[$whichFrags] = $1;
        }
      }
    }
    elsif($inGapSizes && $newline =~ /Number of negatives: (\d+)/)
    {
      $numNegativeGaps = $1;
    }
    elsif($newline =~ /Number of positives: (\d+)/)
    {
      if($inGapSizes)
      {
        $numPositiveGaps = $1;
        $numGaps = $numPositiveGaps + $numNegativeGaps;
      }
      else
      {
        $numContigs = $1;
      }
    }
    elsif($newline =~ /Min: (\d+).(\d+), Mean: (\d+).(\d+), Max: (\d+).(\d+), Stddev: (\d+).(\d+)/)
    {
      if($inGapSizes)
      {
        $minPositiveGap = $1;
        $meanPositiveGap = $3;
        if($4 >= 50)
        {
          $meanPositiveGap = $meanPositiveGap + 1;
        }
        $maxPositiveGap = $5;
        if($6 >= 50)
        {
          $maxPositiveGap = $maxPositiveGap + 1;
        }
      }
      else
      {
        $minContig = $1;
        $meanContig = $3;
        if($4 >= 50)
        {
          $meanContig = $meanContig + 1;
        }
        $maxContig = $5;
        if($6 >= 50)
        {
          $maxContig = $maxContig + 1;
        }
      }
    }
    elsif($newline =~ /Scaffold graph summary:/)
    {
      last;
    }
  }
  elsif($whichBlock == $inUnitig)
  {
#    print "In unitig\n";
    if($newline =~ /(\d+) reads, (\d+) BAC ends, (\d+) external reads, (\d+) external frags/)
    {
      $numReads = $1;
      $numBacEnds = $2;
      $numExtReads = $3;
      $numExtFrags = $4;
      next;
    }
    elsif($newline =~ /Celera reads:/)
    {
      $whichFrags = $inReads;
      next;
    }
    elsif($newline =~ /BAC ends:/)
    {
      $whichFrags = $inBacEnds;
      next;
    }
    elsif($newline =~ /External reads:/)
    {
      $whichFrags = $inExtReads;
      next;
    }
    elsif($newline =~ /External frags:/)
    {
      $whichFrags = $inExtFrags;
      next;
    }
    elsif($newline =~ /Total:/)
    {
      $whichFrags = $inTotal;
      next;
    }
    elsif($newline =~ /Mateless:/)
    {
      $whichFrags = $inMateless;
      next;
    }
    elsif($newline =~ /External mates:/)
    {
      $whichFrags = $inExternal;
      next;
    }

    # detect which line in the block
    if($whichFrags == $inMateless)
    {
      if($newline =~ /(\d+) reads/)
      {
        $matelessUnitig[$inReads] = $1;
      }
      elsif($newline =~ /(\d+) BAC ends/)
      {
        $matelessUnitig[$inBacEnds] = $1;
      }
      elsif($newline =~ /(\d+) external reads/)
      {
        $matelessUnitig[$inExtReads] = $1;
      }
      elsif($newline =~ /(\d+) total/)
      {
        $matelessUnitig[$inTotal] = $1;
      }
    }
    elsif($whichFrags == $inExternal && \
          $newline =~ /(\d+) fragments with external mates/)
    {
      $interUnitig = $1;
    }
    else
    {
      if($newline =~ /(\d+) happy intra-unitig/)
      {
        $happyUnitig[$whichFrags] = $1;
      }
      elsif($newline =~ /(\d+) mis-oriented intra-unitig/)
      {
        $misoUnitig[$whichFrags] = $1;
      }
      elsif($newline =~ /(\d+) mis-separated too close intra-unitig/)
      {
        $closeUnitig[$whichFrags] = $1;
      }
      elsif($newline =~ /(\d+) mis-separated too far intra-unitig/)
      {
        $farUnitig[$whichFrags] = $1;
      }
    }
    next;
  }
  elsif($whichBlock == $inContig)
  {
#    print "In contig\n";
    # detect the block
    if($newline =~ /(\d+) reads, (\d+) BAC ends, (\d+) external reads, (\d+) external frags/)
    {
      $numReads = $1;
      $numBacEnds = $2;
      $numExtReads = $3;
      $numExtFrags = $4;
      next;
    }
    elsif($newline =~ /Celera reads:/)
    {
      $whichFrags = $inReads;
      next;
    }
    elsif($newline =~ /BAC ends:/)
    {
      $whichFrags = $inBacEnds;
      next;
    }
    elsif($newline =~ /External reads:/)
    {
      $whichFrags = $inExtReads;
      next;
    }
    elsif($newline =~ /External frags:/)
    {
      $whichFrags = $inExtFrags;
      next;
    }
    elsif($newline =~ /Total:/)
    {
      $whichFrags = $inTotal;
      next;
    }
    elsif($newline =~ /Mateless:/)
    {
      $whichFrags = $inMateless;
      next;
    }
    elsif($newline =~ /External mates:/)
    {
      $whichFrags = $inExternal;
      next;
    }

    if($newline =~ /Unitig sizes:/)
    {
      $inUnitigSizes = 1;
      next;
    }
    elsif($newline =~ /Surrogate sizes:/)
    {
      $inUnitigSizes = 0;
      next;
    }
    
    if($newline =~ /Number of positives: (\d+)/)
    {
      if($inUnitigSizes)
      {
        $numUnitigs = $1;
      }
      else
      {
        $numSurrogates = $1;
      }
      next;
    }
    elsif($newline =~ /Min: (\d+).(\d+), Mean: (\d+).(\d+), Max: (\d+).(\d+), Stddev: (\d+).(\d+)/)
    {
      if($inUnitigSizes)
      {
        $minUnitig = $1;
        $meanUnitig = $3;
        if($4 >= 50)
        {
          $meanUnitig = $meanUnitig + 1;
        }
        $maxUnitig = $5;
        if($6 >= 50)
        {
          $maxUnitig = $maxUnitig + 1;
        }
      }
      else
      {
        $minSurrogate = $1;
        $meanSurrogate = $3;
        if($4 >= 50)
        {
          $meanSurrogate = $meanSurrogate + 1;
        }
        $maxSurrogate = $5;
        if($6 >= 50)
        {
          $maxSurrogate = $maxSurrogate + 1;
        }
      }
      next;
    }

    # detect which line in the block
    if($whichFrags == $inMateless)
    {
      if($newline =~ /(\d+) reads/)
      {
        $matelessContig[$inReads] = $1;
      }
      elsif($newline =~ /(\d+) BAC ends/)
      {
        $matelessContig[$inBacEnds] = $1;
      }
      elsif($newline =~ /(\d+) external reads/)
      {
        $matelessContig[$inExtReads] = $1;
      }
      elsif($newline =~ /(\d+) total/)
      {
        $matelessContig[$inTotal] = $1;
      }
    }
    elsif($whichFrags == $inExternal && $newline =~ /(\d+) fragments with external mates/)
    {
#      print "Switching from contig to scaffold\n";
      $interContig = $1;
      $whichBlock = $inScaffold;
    }
    else
    {
      if($newline =~ /(\d+) happy intra-contig/)
      {
        $happyContig[$whichFrags] = $1 - $happyUnitig[$whichFrags];
      }
      elsif($newline =~ /(\d+) mis-oriented intra-contig/)
      {
        $misoContig[$whichFrags] = $1 - $misoUnitig[$whichFrags];
      }
      elsif($newline =~ /(\d+) mis-separated too close intra-contig/)
      {
        $closeContig[$whichFrags] = $1 - $closeUnitig[$whichFrags];
      }
      elsif($newline =~ /(\d+) mis-separated too far intra-contig/)
      {
        $farContig[$whichFrags] = $1 - $farUnitig[$whichFrags];
      }
    }
  }
}

# done
CloseFiles();
exit(0);
###########################################################


###########################################################
# close input & output files
sub CloseFiles
{
  close(TXTIN);
  close(UMOUT);
  close(CMOUT);
  close(SMOUT);
  close(AMOUT);
  close(HAPPYOUT);
  close(MISOOUT);
  close(FAROUT);
  close(CLOSEOUT);
  close(CTGOUT);
  close(SCFOUT);
  close(FRGOUT);
  close(READOUT);
  close(BEOUT);
  close(EXTROUT);
  close(EXTFOUT);
  close(TYPEOUT);
}
###########################################################


###########################################################
# subroutine to write to any of the mates files
sub ComputeMatePercentages
{
  $numBad = $numMiso + $numClose + $numFar;
  $numMates = $numHappy + $numBad;
  if($numMates > 0)
  {
    $pctHappy = 100 * $numHappy / $numMates;
    $pctMiso = 100 * $numMiso / $numMates;
    $pctClose = 100 * $numClose / $numMates;
    $pctFar = 100 * $numFar / $numMates;
    $pctBad = 100 * $numBad / $numMates;
  }
  else
  {
    $pctHappy = 0;
    $pctMiso = 0;
    $pctClose = 0;
    $pctFar = 0;
    $pctBad = 0;
    $numMates = 0;
  }
}
###########################################################


###########################################################
sub PrintToUnitigMatesFile
{
  $numHappy = $happyUnitig[$allFrags];
  $numMiso = $misoUnitig[$allFrags];
  $numFar = $farUnitig[$allFrags];
  $numClose = $closeUnitig[$allFrags];
  $numMateless = $matelessUnitig[$allFrags];
  $numInter = $interUnitig;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf UMOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $numBad, $numInter, $numMateless, $pctHappy, $pctMiso, $pctClose, $pctFar, $pctBad;
  }
}
###########################################################


###########################################################
sub PrintToContigMatesFile
{
  $numHappy = $happyContig[$allFrags];
  $numMiso = $misoContig[$allFrags];
  $numFar = $farContig[$allFrags];
  $numClose = $closeContig[$allFrags];
  $numMateless = $matelessContig[$allFrags];
  $numInter = $interContig;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf CMOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $numBad, $numInter, $numMateless, $pctHappy, $pctMiso, $pctClose, $pctFar, $pctBad;
  }
}
###########################################################


###########################################################
sub PrintToScaffoldMatesFile
{
  $numHappy = $happyScaffold[$allFrags];
  $numMiso = $misoScaffold[$allFrags];
  $numFar = $farScaffold[$allFrags];
  $numClose = $closeScaffold[$allFrags];
  $numMateless = $matelessScaffold[$allFrags];
  $numInter = $interScaffold;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf SMOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $numBad, $numInter, $numMateless, $pctHappy, $pctMiso, $pctClose, $pctFar, $pctBad;
  }
}
###########################################################


###########################################################
sub PrintToAllMatesFile
{
  $numHappy = $happyUnitig[$allFrags] + $happyContig[$allFrags] + $happyScaffold[$allFrags];
  $numMiso = $misoUnitig[$allFrags] + $misoContig[$allFrags] + $misoScaffold[$allFrags];
  $numFar = $farUnitig[$allFrags] + $farContig[$allFrags] + $farScaffold[$allFrags];
  $numClose = $closeUnitig[$allFrags] + $closeContig[$allFrags] + $closeScaffold[$allFrags];
  $numMateless = $matelessScaffold[$allFrags];
  $numInter = $interScaffold;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf AMOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $numBad, $numInter, $numMateless, $pctHappy, $pctMiso, $pctClose, $pctFar, $pctBad;
  }
}
###########################################################


###########################################################
sub PrintToHappyMatesFile
{
  my $total = $happyUnitig[$allFrags] + $happyContig[$allFrags] + $happyScaffold[$allFrags];
  if($total == 0)
  {
    $pctUnitig = 0;
    $pctContig = 0;
    $pctScaffold = 0;
  }
  else
  {
    $pctUnitig = 100 * $happyUnitig[$allFrags] / $total;
    $pctContig = 100 * $happyContig[$allFrags] / $total;
    $pctScaffold = 100 * $happyScaffold[$allFrags] / $total;
  }

  printf HAPPYOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $happyUnitig[$allFrags], $happyContig[$allFrags], $happyScaffold[$allFrags], $total, $pctUnitig, $pctContig, $pctScaffold;
}
###########################################################


###########################################################
sub PrintToMisoMatesFile
{
  my $total = $misoUnitig[$allFrags] + $misoContig[$allFrags] + $misoScaffold[$allFrags];
  if($total == 0)
  {
    $pctUnitig = 0;
    $pctContig = 0;
    $pctScaffold = 0;
  }
  else
  {
    $pctUnitig = 100 * $misoUnitig[$allFrags] / $total;
    $pctContig = 100 * $misoContig[$allFrags] / $total;
    $pctScaffold = 100 * $misoScaffold[$allFrags] / $total;
  }

  printf MISOOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $misoUnitig[$allFrags], $misoContig[$allFrags], $misoScaffold[$allFrags], $total, $pctUnitig, $pctContig, $pctScaffold;
}
###########################################################


###########################################################
sub PrintToFarMatesFile
{
  my $total = $farUnitig[$allFrags] + $farContig[$allFrags] + $farScaffold[$allFrags];
  if($total == 0)
  {
    $pctUnitig = 0;
    $pctContig = 0;
    $pctScaffold = 0;
  }
  else
  {
    $pctUnitig = 100 * $farUnitig[$allFrags] / $total;
    $pctContig = 100 * $farContig[$allFrags] / $total;
    $pctScaffold = 100 * $farScaffold[$allFrags] / $total;
  }

  printf FAROUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $farUnitig[$allFrags], $farContig[$allFrags], $farScaffold[$allFrags], $total, $pctUnitig, $pctContig, $pctScaffold;
}
###########################################################



###########################################################
sub PrintToCloseMatesFile
{
  my $total = $misoUnitig[$allFrags] + $misoContig[$allFrags] + $misoScaffold[$allFrags] + $farUnitig[$allFrags] + $farContig[$allFrags] + $farScaffold[$allFrags] + $closeUnitig[$allFrags] + $closeContig[$allFrags] + $closeScaffold[$allFrags];
  my $unitigTotal = $misoUnitig[$allFrags] + $farUnitig[$allFrags] + $closeUnitig[$allFrags];
  my $contigTotal = $misoContig[$allFrags] + $farContig[$allFrags] + $closeContig[$allFrags];
  my $scaffoldTotal = $misoScaffold[$allFrags] + $farScaffold[$allFrags] + $closeScaffold[$allFrags];

  if($total == 0)
  {
    $pctUnitig = 0;
    $pctContig = 0;
    $pctScaffold = 0;
  }
  else
  {
    $pctUnitig = 100 * $unitigTotal / $total;
    $pctContig = 100 * $contigTotal / $total;
    $pctScaffold = 100 * $scaffoldTotal / $total;
  }

  printf BADOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $unitigTotal, $contigTotal, $scaffoldTotal, $total, $pctUnitig, $pctContig, $pctScaffold;
}
###########################################################


###########################################################
sub PrintToBadMatesFile
{
  my $total = $closeUnitig[$allFrags] + $closeContig[$allFrags] + $closeScaffold[$allFrags];
  if($total == 0)
  {
    $pctUnitig = 0;
    $pctContig = 0;
    $pctScaffold = 0;
  }
  else
  {
    $pctUnitig = 100 * $closeUnitig[$allFrags] / $total;
    $pctContig = 100 * $closeContig[$allFrags] / $total;
    $pctScaffold = 100 * $closeScaffold[$allFrags] / $total;
  }

  printf CLOSEOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $closeUnitig[$allFrags], $closeContig[$allFrags], $closeScaffold[$allFrags], $total, $pctUnitig, $pctContig, $pctScaffold;
}
###########################################################



###########################################################
sub PrintToScaffoldFile
{
  printf SCFOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $size, $numContigs, $minContig, $meanContig, $maxContig, $numGaps, $numNegativeGaps, $numPositiveGaps, $minPositiveGap, $meanPositiveGap, $maxPositiveGap;
}
###########################################################


###########################################################
sub PrintToContigFile
{
  printf CTGOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $numUnitigs, $minUnitig, $meanUnitig, $maxUnitig, $numSurrogates, $minSurrogate, $meanSurrogate, $maxSurrogate;
}
###########################################################


###########################################################
sub PrintToFrgFile
{
  printf FRGOUT "%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $numReads, $numBacEnds, $numExtReads, $numExtFrags, $numLocales;
}
###########################################################


###########################################################
sub PrintToReadMatesFile
{
  $numHappy = $happyScaffold[$inReads] + $happyContig[$inReads] + $happyUnitig[$inReads];
  $numMiso = $misoScaffold[$inReads] + $misoContig[$inReads] + $misoUnitig[$inReads];
  $numFar = $farScaffold[$inReads] + $farContig[$inReads] + $farUnitig[$inReads];
  $numClose = $closeScaffold[$inReads] + $closeContig[$inReads] + $closeUnitig[$inReads];
  $numMateless = $matelessScaffold[$inReads] + $matelessContig[$inReads] + $matelessUnitig[$inReads];
  $numInter = $interScaffold;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf READOUT "%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $pctHappy, $pctMiso, $pctClose, $pctFar;
  }
}
###########################################################


###########################################################
sub PrintToBacEndMatesFile
{
  $numHappy = $happyScaffold[$inBacEnds] + $happyContig[$inBacEnds] + $happyUnitig[$inBacEnds];
  $numMiso = $misoScaffold[$inBacEnds] + $misoContig[$inBacEnds] + $misoUnitig[$inBacEnds];
  $numFar = $farScaffold[$inBacEnds] + $farContig[$inBacEnds] + $farUnitig[$inBacEnds];
  $numClose = $closeScaffold[$inBacEnds] + $closeContig[$inBacEnds] + $closeUnitig[$inBacEnds];
  $numMateless = $matelessScaffold[$inBacEnds] + $matelessContig[$inBacEnds] + $matelessUnitig[$inBacEnds];
  $numInter = $interScaffold;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf BEOUT "%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $pctHappy, $pctMiso, $pctClose, $pctFar;
  }
}
###########################################################


###########################################################
sub PrintToExtReadMatesFile
{
  $numHappy = $happyScaffold[$inExtReads] + $happyContig[$inExtReads] + $happyUnitig[$inExtReads];
  $numMiso = $misoScaffold[$inExtReads] + $misoContig[$inExtReads] + $misoUnitig[$inExtReads];
  $numFar = $farScaffold[$inExtReads] + $farContig[$inExtReads] + $farUnitig[$inExtReads];
  $numClose = $closeScaffold[$inExtReads] + $closeContig[$inExtReads] + $closeUnitig[$inExtReads];
  $numMateless = $matelessScaffold[$inExtReads] + $matelessContig[$inExtReads] + $matelessUnitig[$inExtReads];
  $numInter = $interScaffold;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf EXTROUT "%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $pctHappy, $pctMiso, $pctClose, $pctFar;
  }
}
###########################################################


###########################################################
sub PrintToExtFragMatesFile
{
  $numHappy = $happyScaffold[$inExtFrags] + $happyContig[$inExtFrags] + $happyUnitig[$inExtFrags];
  $numMiso = $misoScaffold[$inExtFrags] + $misoContig[$inExtFrags] + $misoUnitig[$inExtFrags];
  $numFar = $farScaffold[$inExtFrags] + $farContig[$inExtFrags] + $farUnitig[$inExtFrags];
  $numClose = $closeScaffold[$inExtFrags] + $closeContig[$inExtFrags] + $closeUnitig[$inExtFrags];
  $numMateless = $matelessScaffold[$inExtFrags] + $matelessContig[$inExtFrags] + $matelessUnitig[$inExtFrags];
  $numInter = $interScaffold;
  
  ComputeMatePercentages();

  if($numMates > 0)
  {
    printf EXTFOUT "%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\t\%d\n", $scaffoldID, $numHappy, $numMiso, $numClose, $numFar, $pctHappy, $pctMiso, $pctClose, $pctFar;
  }
}
###########################################################


###########################################################
sub PrintToMateTypesFile
{
  my @pctH = (0) x ($allFrags);
  my @pctM = (0) x ($allFrags);
  my @pctF = (0) x ($allFrags);
  my @pctC = (0) x ($allFrags);
  my @pctB = (0) x ($allFrags);
  my $i = 0;
  
  $numHappy = $happyScaffold[$allFrags];
  $numMiso = $misoScaffold[$allFrags];
  $numFar = $farScaffold[$allFrags];
  $numClose = $closeScaffold[$allFrags];
  $numBad = $numMiso + $numFar + $numClose;

  if($numHappy > 0)
  {
    for($i = 0; $i <= $allFrags; $i++)
    {
      $pctH[$i] = 100 * $happyScaffold[$i] / $numHappy;
    }
  }
  if($numMiso > 0)
  {
    for($i = 0; $i <= $allFrags; $i++)
    {
      $pctM[$i] = 100 * $misoScaffold[$i] / $numMiso;
    }
  }
  if($numFar > 0)
  {
    for($i = 0; $i <= $allFrags; $i++)
    {
      $pctF[$i] = 100 * $farScaffold[$i] / $numFar;
    }
  }
  if($numClose > 0)
  {
    for($i = 0; $i <= $allFrags; $i++)
    {
      $pctC[$i] = 100 * $closeScaffold[$i] / $numClose;
    }
  }
  if($numBad > 0)
  {
    for($i = 0; $i <= $allFrags; $i++)
    {
      $pctB[$i] = 100 * ($farScaffold[$i] + $closeScaffold[$i] + $misoScaffold[$i]) / $numBad;
    }
  }
  printf TYPEOUT "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", $scaffoldID, $pctH[$inReads], $pctH[$inBacEnds], $pctH[$inExtReads], $pctH[$inExtFrags], $pctM[$inReads], $pctM[$inBacEnds], $pctM[$inExtReads], $pctM[$inExtFrags], $pctC[$inReads], $pctC[$inBacEnds], $pctC[$inExtReads], $pctC[$inExtFrags], $pctF[$inReads], $pctF[$inBacEnds], $pctF[$inExtReads], $pctF[$inExtFrags], $pctB[$inReads], $pctB[$inBacEnds], $pctB[$inExtReads], $pctB[$inExtFrags];
}
###########################################################


###########################################################
sub WriteToAllFiles
{
#  print "Writing to all files\n";
  PrintToFrgFile();
  PrintToContigFile();
  PrintToScaffoldFile();
  PrintToUnitigMatesFile();
  PrintToContigMatesFile();
  PrintToScaffoldMatesFile();
  PrintToAllMatesFile();
  PrintToHappyMatesFile();
  PrintToMisoMatesFile();
  PrintToFarMatesFile();
  PrintToCloseMatesFile();
  PrintToBadMatesFile();
  PrintToReadMatesFile();
  PrintToBacEndMatesFile();
  PrintToExtReadMatesFile();
  PrintToExtFragMatesFile();
  PrintToMateTypesFile();
}
###########################################################
