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
# $Id: sim4db2tampa.pl,v 1.3 2005-12-16 22:13:07 catmandew Exp $
#

use FileHandle;
use Getopt::Long;

# algo:
# open sim4db file
#   keep track of chromosomes/scaffolds
#   read in mappings into %mappings hashtable
#     identify multiply mapped frags & keep best, if user specified
# close sim4db file
#
# open output files:
#   two for each chromosome/scaffold - intra and inter
#   one libs file
# open frag file
#   as DST messages are encountered, write info to libs file
#   as LKG messages are encountered, look up in %mappings
#     write paired mapping info to appropriate file(s)
# close files

sub usage
{
  printf STDERR "Usage: ${0}  -f fragFilename  -s sim4dbFilename  -o outputPrefix  [-m multiMatches]  [-h]\n" .
    "  -f fragFilename    name of .frg-type file\n" .
    "                       which contains DST and LKG messages\n" .
    "  -s sim4dFilename   name of sim4db output file\n" .
    "                       which lists fragment mappings\n" .
    "  -o outputPrefix    prefix of filenames for output\n" .
    "  -m multiMatches    what to do with fragments with\n" .
    "                       multiple mappings\n" .
    "                       omit this flag to list only uniquely\n" .
    "                         mapped fragments\n" .
    "                       c = list match with highest coverage\n" .
    "                       i = list match with highest identity\n" .
    "                       default is to omit multiply mapped\n" .
    "                         fragments.\n";
  exit;
}

my $fragFilename = "";
my $sim4Filename = "";
my $outputPrefix = "";
my $help = 0;

# "", "c", "i"
my $multiMatch = "";
  
GetOptions("f=s", => \$fragFilename,
           "s=s", => \$sim4Filename,
           "o=s", => \$outputPrefix,
           "m:s", => \$multiMatch,
           "h", => \$help) or die "Option error";

usage if($help != 0 ||
         $fragFilename eq "" || $sim4Filename eq "" || $outputPrefix eq "" ||
         ($multiMatch ne "" && $multiMatch ne "c" && $multiMatch ne "i"));

# user specifies:
#  name of frag file
#  name of sim4db file
#  one of the following:
#    exclude multiply matching fragments
#    include match with highest coverage
#    include match with highest identity
#  prefix of output filenames

# output files - space-delimited
# 1. lib file
#    fields: libID, mean length, stddev
# 2. intra-chromosome/scaffold files
#    one per chromosome/scaffold
#    fields: orientation (I,O,N,A), left fragID, right fragID, libID,
#            coord of 5' end of left frag, coord of 5' end of right frag
# 3. inter-chromosome/scaffold files
#    one per chromosome/scaffold
#    fields: this fragID, this chromosome/sequence ID, this 5' coord,
#            this orientation (A_B, B_A),
#            other fragID, other chromosome/sequence ID, other 5' coord,
#            other orientation
#            libID

# hashtable:
# key = fragID
# value = best mapping with count of mappings, array of
#   chromosome/scaffold ID, 5' coordinate, orientation,
#   identity %, coverage %, number of matches
my %mappings;

# hashtable:
# key = chromosome/scaffold ID
# value = whatever, initially
#   later, set to an array of two filehandles per, intra- and inter-
my %intraFHs;
my %interFHs;

# reusable variables
my $fh;
my @hashArray;
my @fields;

# now read in sim4db file
# format is:
# sim4begin
# fragiid[fraglen-0-0] contigiid[0-0] <numidentities-0-percentid-mappingorient-unknown>
# edef=>frag defline
# ddef=>contig defline
# fragbegin-fragend (contigbegin-contigend) <numidentities-0-percentid>
# sim4end

# mappingorient will be either 'forward' or 'complement'; 'unknown' there
# is the strand orientitation, which makes no sense here, so just ignore
# this. The zeros in contiguid[0-0] should always be zero, the other
# zeros might occasionally be non-zero, just ignore. contig positions are
# always relative to the forward contig. If you see more than one
# alignment block per record, the aligner wanted to insert a splice. I
# generally don't trust these (there shouldn't be that many). The iid's
# start at zero, and fragiid's are consistent between assemblies.

# hashArray will hold these fields:
#   chrom/scaffold, 5' coord, orientation, identity, coverage, #mappings

my $recordLine = 0;
my $linesPerMapping = 6;
my $fragID;
my $i;
my $numMappings = 0;
my $numSplicedMappings = 0;
my $numUniqueContiguousMappings = 0;

my $lineCount = 0;
printf STDERR "Reading sim4db file %s\n", $sim4Filename;
$sim4FH = new FileHandle $sim4Filename, "r" or die "Failed to open $sim4Filename for reading";
while(<$sim4FH>)
{
  $recordLine++;
  s/[\n\r\cZ]//g;

  # line indicating what is mapped to what, but not where
  if($recordLine == 2)
  {
    s/[\[\] <>]/-/g;
    @fields = split "-";

    # printf STDOUT "$_\n";

    # chromosome/scaffold
    $hashArray[0] = $fields[5];
    
    # orientation
    $hashArray[2] = ($fields[13] eq "forward" ? 0 : 1);
    
    # identity (percent)
    $hashArray[3] = $fields[12];
    
    # coverage (bases)
    $hashArray[4] = $fields[10];
    
    # number of matches
    $hashArray[5] = 1;

    $numMappings++;
    next;
  }

  # line indicating which fragment
  if($recordLine == 3)
  {
    s/[\[\] <>=]/-/g;
    @fields = split "-";
    
    $fragID = $fields[2];
    next;
  }
  
  # line indicating where
  if($recordLine == 5)
  {
    s/[\(\) <>]/-/g;
    @fields = split "-";
    
    # printf STDOUT "$_\n";
    
    $hashArray[1] = ($hashArray[2] == 0 ? $fields[3] : $fields[4]);
    next;
  }

  # end of record marker
  if($_ eq "sim4end")
  {
    printf(STDERR "\r%10d", $lineCount) if(++$lineCount % 10000 == 0);
  
    # check to omit mappings with multiple spliced alignments
    if($recordLine == $linesPerMapping)
    {
      if(defined($mappings{$fragID}))
      {
        $numUniqueContiguousMappings-- if($mappings{$fragID}[5] == 1);
        if($multiMatch eq "")
        {
          $mappings{$fragID}[5]++;
          next;
        }

        $numMatches = $mappings{$fragID}[5] + 1;

        # if user specified identity preference,
        # replace this match if it has higher identity or
        # equal identity and higher coverage
        $mappings{$fragID} = [@hashArray]
          if($multiMatch eq "i" && ($mappings{$fragID}[3] < $hashArray[3] ||
                                    ($mappings{$fragID}[3] == $hashArray[3] &&
                                     $mappings{$fragID}[4] < $hashArray[4])));
        
        # if user specified coverage preference,
        # replace this match if it has higher coverage or
        # equal coverage and higher identity
        $mappings{$fragID} = [@hashArray]
          if($multiMatch eq "c" && ($mappings{$fragID}[4] < $hashArray[4] ||
                                    ($mappings{$fragID}[4] == $hashArray[4] &&
                                     $mappings{$fragID}[3] < $hashArray[3])));

        # whether switched or not, update the number of matches
        $mappings{$fragID}[5] = $numMatches;
      }
      else
      {
        $numUniqueContiguousMappings++;
        $mappings{$fragID} = [@hashArray];
      }
    }
    else
    {
      $numSplicedMappings++;
    }
    $recordLine = 0;
  }
}
close($sim4FH);
printf(STDERR "\r%10d\n", $lineCount);

# unique, contiguous mappings

printf STDOUT "\nRead %9d total mappings\n", $numMappings;
printf STDOUT "     %9d were spliced\n", $numSplicedMappings;
printf STDOUT "     %9d were unique and contiguous\n", $numUniqueContiguousMappings;

my $inDST = 0;
my $inLKG = 0;

my @dst;
my %dstFields = ("acc" => 0,
                 "mea" => 1,
                 "std" => 2,);

my @lkg;
my %lkgFields = ("fg1" => 0,
                 "fg2" => 1,
                 "dst" => 2,);

my $numLibs = 0;
my $numMatePairs = 0;

# read fragments
printf STDERR "Reading $fragFilename\n";
my $fragsFH = new FileHandle $fragFilename, "r" or
  die "Failed to open $fragFilename for reading";

# write libs
my $fn = $outputPrefix . "Libs.txt";
printf STDERR "Writing libraries to $fn\n";
$libsFH = new FileHandle $fn, "w" or die "Failed to open $fn for writing";

while(<$fragsFH>)
{
  s/[\n\r\cZ]//g;

  printf(STDERR "\r%10d", $i) if(++$i % 10000 == 0);
  
  # {DST line starts a library
  if(index($_, "{DST") == 0)
  {
    $numLibs++;
    $inLKG = 0;
    $inDST = 1;
    next;
  }
  if($inDST == 1)
  {
    # action field - must be add
    if(index($_, "act:") == 0)
    {
      $inDST = 0 if(index($_, "act:A") < 0);
      next;
    }

    # end of record - write to libs file
    if(index($_, "}") == 0)
    {
      printf $libsFH "%s %f %f\n", $dst[0], $dst[1], $dst[2];
      $inDST = 0;
      next;
    }

    # otherwise it's a field to hold onto
    @fields = split ":";
    $dst[$dstFields{$fields[0]}] = $fields[1];
  }

  # {LKG line starts a mate pair
  if(index($_, "{LKG") == 0)
  {
    $numMatePairs++;
    $inLKG = 1;
    $inDST = 0;
    next;
  }
  
  if($inLKG == 1)
  {
    # action field - must be add
    if(index($_, "act:") == 0)
    {
      $inLKG = 0 if(index($_, "act:A") < 0);
      next;
    }

    # ignore certain fields
    next if(index($_, "typ:") == 0 || index($_, "etm:") == 0);

    # end of record - process wrt mappings
    if(index($_, "}") == 0)
    {
      # print the data to the appropriate file(s)

      # don't print if either fragment is unmapped or
      # if one has multiple mappings and user doesn't want them
      next if(!defined($mappings{$lkg[$lkgFields{"fg1"}]}) ||
              !defined($mappings{$lkg[$lkgFields{"fg2"}]}) ||
              ($multiMatch eq "" &&
               ($mappings{$lkg[$lkgFields{"fg1"}]}[5] > 1 ||
                $mappings{$lkg[$lkgFields{"fg2"}]}[5] > 1)));

      # printf STDERR "Not skipping this LKG\n";
      # if both frags are on same chrom, write to intra file
      if($mappings{$lkg[$lkgFields{"fg1"}]}[0] ==
         $mappings{$lkg[$lkgFields{"fg2"}]}[0])
      {
        my $orient;
        my $left;
        my $right;

        # determine which frag is left vs right
        if($mappings{$lkg[$lkgFields{"fg1"}]}[1] <
           $mappings{$lkg[$lkgFields{"fg2"}]}[1])
        {
          $left = "fg1";
          $right = "fg2";
        }
        else
        {
          $left = "fg2";
          $right = "fg1";
        }

        # determine orientation of pair
        if($mappings{$lkg[$lkgFields{$left}]}[2] == 0)
        {
          # left frag is forward
          $orient = ($mappings{$lkg[$lkgFields{$right}]}[2] == 0 ? "N" : "I");
        }
        else
        {
          # left frag is complemented
          $orient = ($mappings{$lkg[$lkgFields{$right}]}[2] == 0 ? "O" : "A");
        }

        # intra:
        #   orientation (I,O,N,A), leftID, rightID, libID, left 5', right 5'
        my $mode = 
            (defined($intraFHs{$mappings{$lkg[$lkgFields{"fg1"}]}[0]}) ?
             "a" : "w");
        $intraFHs{$mappings{$lkg[$lkgFields{"fg1"}]}[0]} = 1;
        $fn = $outputPrefix . "_" .
          $mappings{$lkg[$lkgFields{"fg1"}]}[0] .
          "_intra.txt";
        my $thisFH =
            new FileHandle $fn, $mode or die "Failed to open $fn for writing";
        printf($thisFH "%s %s %s %s %u %u\n",
               $orient,
               $lkg[$lkgFields{$left}],
               $lkg[$lkgFields{$right}],
               $lkg[$lkgFields{"dst"}],
               $mappings{$lkg[$lkgFields{$left}]}[1],
               $mappings{$lkg[$lkgFields{$right}]}[1]);
        close($thisFH);
      }
      else
      {
        # write to two inter-chrom files
        
        # inter - first set of info is for frag on 'this' chrom:
        #   fragID, chromID, frag 5', orientation (A_B, B_A),
        #   otherFragID, chromID, frag 5', orientation, libID
        my $fi;
        for($fi = 0; $fi < 2; $fi++)
        {
          my $mode = 
            (defined($interFHs{$mappings{$lkg[$fi]}[0]}) ? "a" : "w");
          $interFHs{$mappings{$lkg[$fi]}[0]} = 1;
          $fn = $outputPrefix . "_" .
            $mappings{$lkg[$fi]}[0] .
            "_inter.txt";
          my $thisFH = 
            new FileHandle $fn, $mode or die "Failed to open $fn for writing";
          printf($thisFH "%s %s %u %s %s %s %u %s %s\n",
                 $lkg[$fi],
                 $mappings{$lkg[$fi]}[0],
                 $mappings{$lkg[$fi]}[1],
                 ($mappings{$lkg[$fi]}[2] == 0 ? "A_B" : "B_A"),
                 $lkg[1-$fi],
                 $mappings{$lkg[1-$fi]}[0],
                 $mappings{$lkg[1-$fi]}[1],
                 ($mappings{$lkg[1-$fi]}[2] == 0 ? "A_B" : "B_A"),
                 $lkg[2]);
          close($thisFH);
        }
      }
      $inLKG = 0;
      next;
    }

    # otherwise it's a field to hold onto
    @fields = split ":";
    # skip if not innie
    if($fields[0] eq "ori")
    {
      $inLKG = 0 if($fields[1] ne "I");
      next;
    }
    $lkg[$lkgFields{$fields[0]}] = $fields[1];
  }
}
close($fragsFH);
printf(STDERR "\r%10d\n", $i);
printf(STDOUT "\n%10d matepairs in frg file\n", $numMatePairs);

close($libsFH);
