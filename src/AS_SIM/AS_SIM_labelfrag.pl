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
#ASM_SIM_labelfrag.pl
#   This script is invoked by AS_SIM_celsim.awk.
#   Inputs:
#       <cmsfilename>     Comment file from celsim invocation
#       <frgfilename>     Raw fragment file (either a urc or frg file)
#       .note.<Seed>       File listing elements for annotation
#   Outputs:
#       <frgfilename>     Modified in-place with repeat annotations
#
#   Invocation:
#       perl ASM_SIM_labelfrag.pl <cmsfilename> <frgfilename> <noteFileName>
#
#  Annotations are output in the src field of relevant fragments,
#  and or of the format:
#    <Elem>.<Regen_Nat> [<asp_Nat>,<aep_Nat>] [<fragsp_Nat>,<fragep_Nat>] 
#
#  <Elem> is the letter specifying a DNA sequence element, 
#  <Regen_nat> is the index of regeneration of the DNA sequence element,
#  asp_Nat and aep_Nat are the range of coordinates of the <Elem> bases 
#       that appear in this fragment in the coordinates of the DNA 
#       sequence element (5-prime to 3-prime).
#  fragsp_Nat and fragep_Nat are the range of the fragment coordinates 
#       (5-prime to 3-prime) that overlaps the given element instance.
#
#  Note that in the output if the <fragsp_Nat> > <fragep_Nat>, then the 
#  DNA is Watson-Crick complemented.
#  The input files uses [rf] to specify the orientation of the DNA
#  instance in the global genome coordinates.
#
# Operation:
#   The .cms file is read, an array of interesting instances is constructed, 
#   ordered by the left endpoint.  Binary search is used to find the range of 
#   interesting instances that may intersect the fragment. Note that it is
#   assumed that the instances are already in sorted order in the cms file.
#
#  Changed by Knut Reinert 06/04/99
# Added to the source field some information :
# Instead of
#  <Elem>.<Regen_Nat> [<asp_Nat>,<aep_Nat>] [<fragsp_Nat>,<fragep_Nat>] 
# the output is now
#  <Elem>.<Regen_Nat>.<occurence> (noOfOccurences) [<asp_Nat>,<aep_Nat>] [<fragsp_Nat>,<fragep_Nat>] 
# where occurence means the number in the left to right order in Celsim
# coordinates. Behind that the total number of occurences of this element is given.


$verbose = 0;
$comfilename = shift(@ARGV);
$frgfilename = shift(@ARGV);
$notfile = shift(@ARGV);

# print STDERR "comFilename = $comfilename   frgFileName = $frgfilename notes = $notfile\n";


# Read the annotation file
#
if(readNot($notfile) == 0){   # no annotations -- we're done
    print STDERR "# No Annotations -- leaving labelfrag\n";
    exit;
}

readComs($comfilename ); # Read the .cms file
EditFrags($frgfilename ); # Edit the .frg file



#################################################################


sub readNot{
    my($notFileName) = @_;


    # If we can't open the notfile, we are done
    if(NULL == open(NOTS, $notFileName )){
	return 0;
    }

    $numNotations = 0;
    %notations = ();  # An associative array
    while(<NOTS>){
	if(/^(([A-Z]))/){
	    $numNotations++;
	    print "\tAnnotating Element $1 \n";
	    if($notations{$1} != NULL){
		die( "***WIERDNESS!!!! $notations{$1}\n");
	    }
	    $notations{$1} = 1;  # Mark this element for notation
	}
    }
    return $numNotations;
}


#
#EditFrags(fileName)
#
#   Uses in-place editing to transform annotate each FRG in the .frg file
#   With information on repeat content.
#
sub EditFrags{
    my($fragFileName) = @_;
    my $filenamebak, $cnt, $i, $tf, $esp, $eep, $asp, $aep, $left, $right, $startSrc;

    @ARGV = ($fragFileName);
    $^I = "";  # Nuke the old file after the edit
# first test that we can open for r/w the file
    open(FRAG, "$fragFileName") ||
	die "Sorry...cannot open $fragFileName";
    close(FRAG);
#    print STDERR "EditFrags... $fragFileName $filenamebak\n";

    $cnt = 0;
    $bStart = -1;
# Read each line of the file
    while(<>){
	$cnt++;
	if($startSrc == 0){  #Print each line, note the presence of a src:
	    if(/^src:/){
		$startSrc = 1;
	    }
	    print ;
	}elsif($startSrc == 1){ # Inside a src field, do special stuff

	    if(/^\[([0-9]*),([0-9]*)/){ # find the interval of the fragment.
		$bStart = $1; # 
		$bEnd = $2;
#			print STDERR "bStart = $bStart bEnd = $bEnd\n";
		print ;
		
	    } elsif(/^\./){ # find the end of the field
		$startSrc = 0;
                if($bStart > 0){
		$sp = min($bStart,$bEnd);
		$ep = max($bStart,$bEnd);
  	        # print STDERR "findInstances $sp, $ep\n";
		# Use a binary search to find the instances that overlap
		# this fragment.
		($left, $right) = findInstances($sp,$ep, 1,$instances ); 

		# Output the guilty elements 
		#
		if($left > 0 && $right >0){
#		    print STDERR "After findInstance $left, $right\n";
		    while($left <= $right){
			($tf, $esp, $eep, $asp, $aep) = 
			    intervalsOverlap(@start[$left], @end[$left],
					     @direction[$left],
					     $bStart, $bEnd);
#			print STDERR  "$element[$left]\.$regen[$left] [$esp,$eep] [$asp,$aep] \n";
			print  "$element[$left]\.$regen[$left]\.$invocation{$start[$left].$end[$left]} ($count{$element[$left].$regen[$left]}) [$esp,$eep] [$asp,$aep] \n";
			$left++;
		    }
		  }
		$bStart = -1;
	      }
		print  ".\n";
	    } else{
		print ;
	    }

	}
    }
}    


#
# readComs(filename)
#
# read the .cms file, and build the tree of instances
#
#

sub readComs{
    my($comsFileName) = @_;

    print STDERR "Opening $comsFileName\n";
    open(COMS, "$comsFileName") ||
	die "Sorry...cannot open $comsFileName";
    $Root = "0";
#    print STDERR "Root = $Root \n";
    $instances = 0;
    while(<COMS>){
#	print STDERR "Read $_\n";
      if(/^#  ([A-Z])\.0\s*\(([0-9]*)\)/) {  # The Root instance (last defined)
	 if($verbose){
 	   print STDERR "Root Element  $1 \n";
         }
	 if($notations{$1}){
	  $instances++;
	  @level[$instances] = 0;
	  @type[$instances] = " "; # basis/global/concat
	  @element[$instances] = $1; # letter
	  @regen[$instances] = 0;  # regen index
	  @direction[$instances] = "f"; # r or f
	  @start[$instances] = 0;     # interval
	  @end[$instances] = 2;   
          #      @ElementLines[$instances++] = $_;
      }

     }
      # Instances within the file
      if(/^#(\s*)([>=]) ([A-Z])\.([0-9]*)([rf]) at ([0-9]*)-([0-9]*)/){
	 if($verbose){
	     print STDERR "Element $3 \n";
	 }
	 if($notations{$3}){
	     $instances++;
	     @level[$instances] = length($1);
	     @type[$instances] = $2; # basis/global/concat
	     @element[$instances] = $3; # letter
	     @regen[$instances] = $4;  # regen index
	     @direction[$instances] = $5; # r or f
	     @start[$instances] = $6;     # interval
	     @end[$instances]   = $7;      
	     if( $count{$3.$4} == undef ){
	       $count{$3.$4} = 1;
	       $invocation{$6.$7} = 1;
	     }
	     else{
	       $count{$3.$4}++;
	       $invocation{$6.$7} = $count{$3.$4};
	     }


	     
#	      print STDERR "param 5 = $5";
#	  if($5 eq "f") {
#	      print STDERR " forward\n";
#	  } else {
#	      if ($5 eq "r") {
#		  print STDERR " reversed\n";
#		  @start[$instances] = $7;     # interval
#		  @end[$instances]   = $6;       
#	      } else {
#		  die("A fragment that is neither forward or reversed.\n");
#	      }
#	  }
	  if($verbose){
	      print STDERR "mychildren = @$mychildren\n";
	  }
#      @ElementLines[$instances] = $_;
	 }
     }
  }
# Now print the digested version of the lines
#
    if($verbose == 1){
	    print STDERR " instances = $instances \n";
	for($i = 0 ; $i <= $instances; $i++){
	    if( @direction[$i] ne "r" ) {
		print STDERR " $i @level[$i] @element[$i].@regen[$i] (@start[$i],@end[$i]) \n";
	    } else {
		print STDERR " $i @level[$i] @element[$i].@regen[$i] (@end[$i],@start[$i]) \n";
	    } 
	}
    }
}


#
# (true/false, aosp, aoep, bosp, boef) = intervalsOverlap(asp,aep,adr,bsp,bep)
#
#  Compute whether there is an overlap, and the intervals of the two pieces
#  that overlap. The outputs are relative coordinates of the overlap.
#
sub intervalsOverlap{
    my ($aStart, $aEnd, $aDir, $bStart, $bEnd) = @_;
    my $aep, $asp, $bep, $bsp, $tf, $tmp, $aFragmentReversed, $bFragmentReversed;

#    print STDERR "internalsOverlap ($aStart, $aEnd) ($bStartp, $bEndp) ";
#   From the parameters, find the first normalize intervals
    $aFragmentReversed = 0;
    $bFragmentReversed = 0;

    if($aDir eq "r") {
	$aFragmentReversed = 1 - $aFragmentReversed;
    }

    if($aStart > $aEnd){
	$aFragmentReversed = 1;
	$tmp = $aStart;
	$aStart = $aEnd;
	$aEnd = $tmp;
    }
    if($bStart > $bEnd){
	$bFragmentReversed = 1;
	$tmp = $bStart;
	$bStart = $bEnd;
	$bEnd = $tmp;
    }

    if($aEnd < $bStart ||
       $bEnd < $aStart){
	$tf = 0;
	return($tf,0,0,0,0);
# No overlap
    }
    

    $tf = 1;
    if( ( $aStart <= $bStart && $aEnd >= $bEnd)){
	# b contained in a
	$asp = $bStart - $aStart;
	$aep = $bEnd - $aStart;
	$bsp = 0;
	$bep = $bEnd - $bStart;
    }elsif( $bStart <= $aStart && $bEnd >= $aEnd){
	# a contained in b
	$asp = 0;
	$aep = $aEnd - $aStart;
	$bsp = $aStart - $bStart;
	$bep = $aEnd - $bStart;
    }elsif ( $bStart <= $aStart && $bEnd >= $aStart){
	# End of b overlaps start of a
	$asp = 0;
	$aep = $bEnd - $aStart;
	$bsp = $aStart - $bStart;
	$bep = $bEnd - $bStart;
    }elsif ( $aStart <= $bStart && $aEnd >= $bStart){
	# End of a overlaps start of b
	$asp = $bStart - $aStart;  
	$aep = $aEnd - $aStart; # end of a -- suffix
	$bsp = 0;
	$bep = $aEnd - $bStart;
    }else{
#	print "No\n";
	$tf = 0;
    }

    if($aFragmentReversed == 1) {
#       The local coordinates of the fragment are also reversed!
	$asp = ($aEnd - $aStart) - $asp;
	$aep = ($aEnd - $aStart) - $aep;
    }
    if($bFragmentReversed == 1) {
#       The local coordinates of the fragment are also reversed!
	$bsp = ($bEnd - $bStart) - $bsp;
	$bep = ($bEnd - $bStart) - $bep;
    }
    return ($tf, $asp, $aep, $bsp, $bep);
}



# Binary search to find instances that overlap an interval
# Returns a range (first Instance, Last instance) of instances.
# If none found, returns (-1,-1);
#
sub findInstances{
    my($fragStart,$fragEnd, $left, $right) = @_;
    my($midpoint, $leftMost, $rightMost);

    if($left > $right){
	return (-1,-1);
    }
    $midpoint = $left + ceil(($right - $left),2);
#    print STDERR "findInstances $fragStart ($left,$right) midpoint = $midpoint ";

#    print STDERR "(@start[$midpoint], @end[$midpoint])\n";
    if($fragEnd < @start[$midpoint]){
	return findInstances($fragStart,$fragEnd,$left, $midpoint -1);
    }elsif ($fragStart > @end[$midpoint]){
	return findInstances($fragStart,$fragEnd, $midpoint + 1, $right);
    }

#    print STDERR "Scanning for left/rightMost from $midpoint\n";
    # If we get here, then we know that the two intervals intersect
    # We need to scan left, to find the leftmost instance that intersects
    # and scan right to find the rightmost instance that intersects
    $leftMost = $midpoint;
    while($leftMost > 1 &&
	   $fragStart <= @end[$leftMost - 1]){
	# print STDERR "$fragStart >= @end[$leftMost - 1]\n";
	$leftMost--;
    }

    $rightMost = $midpoint;
    while($rightMost < $instances &&
	    $fragEnd >= @start[$rightMost + 1]){
	#  print STDERR "$fragEnd <= @start[rightMost +1]\n";
	$rightMost++;
    }

#    print STDERR "findInstances found leftMost = $leftMost   rightMost = $rightMost\n";
    return ($leftMost, $rightMost);
}




###
### Misc. utility subroutines
###
sub min{
    my($a,$b) = @_;
    if($a <= $b){
	return $a;
    }else{
	return $b;
    }
}

sub max{
    my($a,$b) = @_;
    if($a >= $b){
	return $a;
    }else{
	return $b;
    }
}

sub ceil{
    my($a,$b) = @_;
    my($c);
    if($b == 0){
	# print STDERR "ceil called with $a, $b\n";
	die("ceil called with b == 0\n");
    }
    $c = int $a/$b;
    # print "ceil($a,$b) first result $c\n";
    if($b * $c < $a){
	$c++;
    }
    return $c;
}
    
