#!/usr/local/bin/perl -w
# $Id: reclip.pl,v 1.1 2007-07-03 19:53:59 moweis Exp $
# reclip - Update read vector, quality, and valid clip points based on 
#          results of assembly. 
#
# Written by Martin Shumway
#
 
use strict;
use warnings;
use File::Basename;

# TIGR Modules
use TIGR::Foundation;
use TIGR::AsmLib;

my $tf = new TIGR::Foundation;
my $PRG = $tf->getProgramInfo('name');
my $VERSION="1.10";
my @DEPENDS=("TIGR::Foundation","TIGR::AsmLib");

my $HELPTEXT = qq~
Update read clip points that might have changed as a result of genome assembly.

  reclip    <pull files>  <assembly or retrim files>   -o outprefix
  
  Pull files: 
    prefix.seq.features - Input original read clip points produced by pullfrag.

  Clip files:
    prefix.<modes>      - Input new read clip points where mode can be combo of
      .asm     Assembly file produced by carun
      .tasm    Assembly file produced by run_TA, ca2ta, modContig, slice2tasm, 
                 or AMOS
      .clr     Text file with CLR: seq_name end5 end3
      .clv     Text file with CLV: seq_name end5 end3
      .clb     Text file with CLB: seq_name end5 end3
 
  Options:
    -o <outprefix>  Prefix of the output file, or stdout if not specified. 
    --ref {rid,acc | seq_name,src}  Refer to the read clip point records by id 
          (customary for .asm files), or by alias (customary for .tasm files).  
    --list Produce list files of the following categories:
             outprefix.clr.changed - List of reads whose CLR features changed
             outprefix.clv.changed - List of reads whose CLV features changed
             outprefix.reclip.invalidList of reads whose CLV features changed
    -l <lid> Library EUID identifying which sequencing library [default ALL].
    --metrics  Produce a .metrics file with assembly telemetry related to read
             clip points.  Includes ed_length before and after assembly.
    -N <file>  A file of a subset of read names or rids (see --ref) to reclip

Genome assembly with the Celera Assembler and other tools can produce updates
to read clip points that differ from factory settings.  The reclip program
reconciles these into a single output .bcp file suitable for upload to the
database with aloader.  

Input files are normally supplied by pullfrag as well as the assembly output,
or explicit text files.  Normally all reads mentioned in the .seq.features
file are treated as the considered set, and all those mentioned in the
assembly file are treated as the assembled set.  These sets may be masked off
with the -N option.  This might be used with -l option to perform library
specific analysis of clear range changes. 

The output .bcp file contains a dump of new read clip points for CLR, CLV, CLB 
intervals (clear range, vector free range, and good quality range, respect-
ively).  All coordinates are 1-based inclusive (base-based).  Additional
output files include list files for reclipped reads, a .metrics file suitable
for upload with ametrics, and a outprefix.warnings file that compiles data
exceptions.  The latter should always be inspected upon completion of reclip.

Clip point reconciliation is as follows: If the read's new CLR is not 
contained within the existing CLV, update the CLV end5, end3 or both so that
it does so.  CLB is unaffected.  Do not extend any range outside the length
of the read.  Only those reads mentioned in the assembled set are emitted
in .bcp form for update.

SEE ALSO:  pullfrag, subfrag, retrim, aloader, seqUpload, carun, ametrics
~;

my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;

# =============================== Constants ================================

# Input/Output parameters
my @INPUT_SUFFIXES = ("seq.features", "tasm", "asm");
my %SUPPORTED_INPUTS = ();
map { $SUPPORTED_INPUTS{$_} = 1; } @INPUT_SUFFIXES;
my $REF_BY_RID = 1;
my $REF_BY_ALIAS = 2;
my %REF_MODE = 
(
  "rid"      => $REF_BY_RID,
  "acc"      => $REF_BY_RID,
  "alias"    => $REF_BY_ALIAS,
  "seq_name" => $REF_BY_ALIAS,
  "src"      => $REF_BY_ALIAS,
);

# Operational parameters
my $debug = 0;
my $LOW = 1;
my $MEDIUM = 2;
my $HIGH = 3;
my $warningsfile = "tmp.warnings";
my $AWKFILE = "getClipsFromTasm.awk";
my $AWKSCRIPT = qq~
{
  if (\$1 == "seq_name")
  {
    seq_name = \$2;
  }
  else if (\$1 == "seq_lend")
  {
    seq_lend = \$2;
  }
  else if (\$1 == "seq_rend")
  {
    seq_rend = \$2;
    if (seq_lend > seq_rend)
    {
      print seq_name" "seq_rend" "seq_lend;
    }
    else
    {
      print seq_name" "seq_lend" "seq_rend;
    }
  }
}
~;

# Dispositions, in order of degree of change
my $CONSIDERED = 1;
my $ASSEMBLED  = 2;
my $CHANGED    = 3;
my $RETRIMMED  = 4;
my $RECLIPPED  = 5;

# Trimming parameters
my $PROGRESS_FACTOR = 10000;
my $MINLENGTH = 64;       # Min length of CA read
my $MAXLENGTH = 2048;     # max length of CA read
my $FRG_TOO_SHORT = 1;
my $FRG_TOO_LONG = 2;
my $FRG_OUT_RANGE = 3;
my %FRG_VALID_MESSAGE = 
   (
       $FRG_TOO_SHORT => "less than $MINLENGTH bases.",
       $FRG_TOO_LONG  => "more than $MAXLENGTH bases.",
       $FRG_OUT_RANGE => "clear range outside fragment length.",
   );

# Libraries
my $ALL = "all";


# ========================= Procedures ========================================

# Toss messages to stderr and also to the logfile.  Cannot use 
# regular logfile because of pipeout of this utility.
#
sub progress($)
{
  my ($msg) = @_;
  print STDERR "$msg\n";
  $tf->logLocal($msg, 0, 1);
}

# Obtain new CLR ranges as computed by the CA assembly 
# INputs:
#   infile - Input .asm filename
#   rh_NewClr  - Ref to output lookup table of new clear ranges
#   rh_Emit - Ref to considered set of reads
#   rh_Rid2Alias - Read EUID to Alias lookup
#
sub getClipsFromAsm($$$$)
{
  my ($infile, $rh_NewClr, $rh_Emit, $rh_Rid2Alias) = @_;

  my $if = new IO::File("< $infile") or $tf->bail("Cannot open input file \'$infile\' ($!)");
   
  my $nassembled = 0; 
  while (my $rec = getCARecord($if))
  {
    my ($type, $rh_fields, $recs) = parseCARecord($rec);
    my %fields = %{$rh_fields}; 
    if ($type eq "AFG")
    {
      my $accS = $fields{acc};
      $accS =~ m/\((\d+),\d+\)/;
      my $acc = $1;

      # Presume that the read will either exist in the considered set (Emit table) as
      # a native EUID (rid), or as a seq_name (alias).  Either hit means it's in the
      # considered set and will be extracted from the assembly.
      #
      my $rid = $acc; 
      my $alias = ${$rh_Rid2Alias}{$rid};
      next if ( !(exists ${$rh_Emit}{$rid}  ||  exists ${$rh_Emit}{$alias}) );

      my $clrS = $fields{clr};
      my ($clrl, $clrr) = split /,/,$clrS;
      ${$rh_NewClr}{$rid} = [$clrl + 1, $clrr];  # Convert to 1-based inclusive (base-based)
      ++$nassembled;
    }
  }

  progress("Phase 2: $nassembled records scanned.") if ($nassembled > 1 && $nassembled % $PROGRESS_FACTOR == 1); 
  $if->close();

  return $nassembled;
}

sub getClipsFromTasm($$$$)
{
  my ($infile, $rh_NewClr, $rh_Emit, $rh_Alias2Rid) = @_;

  my $tf = new IO::File("> $AWKFILE") or $tf->bail("Failed to open .awk file \'$AWKFILE\' ($!))");
  $tf->print($AWKSCRIPT);
  $tf->close() or $tf->bail("Failed to create .awk file \'$AWKFILE\' ($!))");

  my $TASMPIPE = "cat $infile | grep seq_ | grep -v seq_id | awk -f $AWKFILE | sort |"; 
  my $if = new IO::File($TASMPIPE) or $tf->bail("Cannot open input pipe cmd \'$TASMPIPE\' ($!)");
  progress("Phase 1 : Scanning input file \'$infile\'...");
   
  my $nassembled = 0; 
  while (my $line = $if->getline())
  {
    chop $line;
    my @f = split /\s+/,$line;

    # Presume that the read will either exist in the considered set (Emit table) as
    # a native EUID (rid), or as a seq_name (alias).  Either hit means it's in the
    # considered set and will be extracted from the assembly.
    #
    my $alias = $f[0];
    my $rid = ${$rh_Alias2Rid}{$alias}; 
    next if ( !(exists ${$rh_Emit}{$rid}  ||  exists ${$rh_Emit}{$alias}) );
    
    ${$rh_NewClr}{$rid} = ($f[1] < $f[2])? [$f[1], $f[2]] : [$f[2], $f[1]];  # Assume 1-based inclusive (base-based)

    ++$nassembled;
  }
 
  unlink $AWKFILE unless ($debug > 0);
  progress("Phase 2: $nassembled records scanned.") if ($nassembled > 1 && $nassembled % $PROGRESS_FACTOR == 1); 

  return $nassembled;
}

sub getClipsFromCsv($$$$)
{
  my ($infile, $rh_Table, $rh_Emit, $rh_Rid2Alias) = @_;

  my $if = new IO::File("< $infile") or $tf->bail("Cannot open input file \'$infile\' ($!)");

  my $nrecords = 0;
  while (my $line = $if->getline())
  {
    chop $line;
    my @f = split /\s+/,$line;
    ${$rh_Table}{$f[0]} = [$f[1], $f[2]];  # Assume 1-based inclusive (base-based)

    ++$nrecords;
  }

  $if->close();

  return $nrecords;
}

sub computeNewCLV($$$$)
{
  my ($new_clrl, $new_clrr, $old_clvl, $old_clvr) = @_;

  my $disposition = $CHANGED;
  my $new_clvl = undef;
  my $new_clvr = undef;

  if    ($new_clrl < $old_clvl && $new_clrr > $old_clvr) 
  {
    # New contains the old, reset both clvl and clvr
    $new_clvl = $new_clrl;
    $new_clvr = $new_clrr;
    $disposition = $RECLIPPED;
  }
  elsif ($new_clrr < $old_clvl) 
  {
    # left exclusion, reset clvl and clvr
    $new_clvl = $new_clrl;
    $new_clvr = $new_clrr;
    $disposition = $RECLIPPED;
  }
  elsif ($new_clrl < $old_clvl && $new_clrr <= $old_clvr) 
  {
    # left intersect, extend clvl
    $new_clvl = $new_clrl;
    $new_clvr = $old_clvr;
    $disposition = $RECLIPPED;
  }
  elsif ($new_clrl >= $old_clvl && $new_clrr <= $old_clvr)
  {
    # contained, leave old values
    $new_clvl = $old_clvl;
    $new_clvr = $old_clvr;
    $disposition = $RETRIMMED;
  }
  elsif ($new_clrl <= $old_clvr && $new_clrr > $old_clvr)
  {
    # right intersect, extend clvr
    $new_clvl = $old_clvl;
    $new_clvr = $new_clrr;
    $disposition = $RECLIPPED;
  }
  elsif ($new_clrl > $old_clvr)
  {
    # right exclusion, reset clvl and clvr
    $new_clvl = $new_clrl;
    $new_clvr = $new_clrr;
    $disposition = $RECLIPPED;
  }
  else
  {
    # should never get here
    $tf->bail("Assertion: Unsupported coordinates.");
  }

  return ($new_clvl, $new_clvr, $disposition);
}

sub toBCP($$$$$$)
{
  my ($outfilename, $rh_Emit, $rh_NewClb, $rh_NewClv, $rh_NewClr, $rh_Rid2Alias) = @_;

  my $of = new IO::File("> $outfilename") or $tf->bail("Cannot open .bcp file \'$outfilename\' ($!)"); 

  foreach my $rid (sort keys %{$rh_Emit})
  {
    next if (${$rh_Emit}{$rid} < $CHANGED);

    # Compile the fields of the extended BCP format
    my $seq_name = (exists ${$rh_Rid2Alias}{$rid})? ${$rh_Rid2Alias}{$rid} : $rid;
    my $new_clbl = ${$rh_NewClb}{$rid}->[0];
    my $new_clbr = ${$rh_NewClb}{$rid}->[1];
    my $new_clnl = ""; 
    my $new_clnr = "";
    my $new_clrl = ${$rh_NewClr}{$rid}->[0];
    my $new_clrr = ${$rh_NewClr}{$rid}->[1];
    my $new_clvl = ${$rh_NewClv}{$rid}->[0];
    my $new_clvr = ${$rh_NewClv}{$rid}->[1];
    my $new_clzl = ""; 
    my $new_clzr = "";
    my $new_seq = "";
    my $new_qual = "";
    my $new_pos = "";
    my $mod_date = $tf->getSybaseDate();   
    $mod_date =~ s/'//g;
    my $id = 0;
    my $line = $seq_name . "\t" .
               $new_clbl . "\t" .
               $new_clbr . "\t" .
               $new_clnl . "\t" .
               $new_clnr . "\t" .
               $new_clrl . "\t" .
               $new_clrr . "\t" .
               $new_clvl . "\t" .
               $new_clvr . "\t" .
               $new_clzl . "\t" .
               $new_clzr . "\t" .
               $new_seq  . "\t" .
               $new_qual . "\t" .
               $new_pos  . "\t" .
               $mod_date . "\t" .
               $id       . "\n" 
               ; 
    $of->print($line) or $tf->bail("Failed to write output .bcp file \'$outfilename\' ($!)");
  } 
  $of->close() or $tf->bail("Failed to close output .bcp file \'$outfilename\' ($!)");
}

sub getClipsFromSeqFeatures($$$$$$$$$)
{
  my ($infile, $rh_OldClv, $rh_OldClb, $rh_OldClr, $filter, $rh_FilterSet, $rh_Rid2Alias, $rh_Alias2Rid, $rh_Emit) = @_;

  my $if = new IO::File("< $infile") or $tf->bail("Cannot open input file \'$infile\' ($!)");

  my $nconsidered = 0;
  while (my $line = $if->getline())
  {
    chop $line;
    my @f = split /\s+/,$line;

    # If there is a filter subset, then only obtain data from the .seq.features file
    # that pertains to the subset (the considered set).
    my $rid = $f[0];
    my $alias = $f[1];
    ${$rh_Rid2Alias}{$f[0]} = $f[1];  # Record the rid to Alias lookup for those reads that will be munged
    ${$rh_Alias2Rid}{$f[1]} = $f[0];  # Record the Alias to Rid lookup for those reads that will be munged
    next if ($filter  &&  ( !(exists ${$rh_FilterSet}{$rid}  ||  exists ${$rh_FilterSet}{$alias}) ) );
    
    ${$rh_OldClv}{$rid} = [$f[4], $f[5]];  # Assume 1-based inclusive (base-based)
    ${$rh_OldClb}{$rid} = [$f[6], $f[7]];  # Assume 1-based inclusive (base-based)
    ${$rh_OldClr}{$rid} = [$f[2], $f[3]];  # Assume 1-based inclusive (base-based)

    ${$rh_Emit}{$rid} = $CONSIDERED;
    ++$nconsidered;
  }

  $if->close();

  return $nconsidered;
}

# toMetrics
#
# Print the following report:
#  [Reclip]
#    TotalReadsConsidered
#    TotalReadsReclipped
#    TotalReadsExtended
#    TotalReadsContracted
#    TotalReadsUnchanged
#    TotalReadsExchanged
#    TotalReadsExtended5PrimeOnly
#    TotalReadsExtended3PrimeOnly
#    TotalReadsExtended5Prime3Prime
#    TotalReadsContracted5PrimeOnly
#    TotalReadsContracted3PrimeOnly
#    TotalReadsContracted5Prime3Prime
#    StartingMeanLength
#    StartingMedianLength
#    StartingMinLength
#    StartingMaxLength
#    StartingStdevLength
#    EndingMeanLength
#    EndingMedianLength
#    EndingMinLength
#    EndingMaxLength
#    EndingStdevLength
#
#
sub toMetrics($$$$$$$$$$)
{
  my ($outfilename, $rh_Emit, $rh_OldClv, $rh_OldClr, $rh_NewClv, $rh_NewClr, $lid, $nconsidered, $nassembled, $nchanged) = @_;

  # Compute aggregate assembly totals
  my $totalReadsConsidered = $nconsidered; 
  my $totalReadsAssembled  = $nassembled;
  my $totalReadsChanged    = $nchanged;

  # Sweep through assembled reads and compute running totals.
  my $aggregateClearRangeConsidered = 0;
  my $aggregateVectorFreeConsidered = 0;
  my $aggregateClearRangeAssembled  = 0;
  my $aggregateVectorFreeAssembled  = 0;
  my $aggregateChange               = 0;                               
  my $aggregateExpansion            = 0;                               
  my $aggregateContraction          = 0;                               
  my $aggregate5PrimeExpansion      = 0;                               
  my $aggregate3PrimeExpansion      = 0;                               
  my $aggregate5PrimeContraction    = 0;                               
  my $aggregate3PrimeContraction    = 0;                               
  my $total5PrimeExpanded           = 0;                     
  my $total5PrimeContracted         = 0;                     
  my $total3PrimeExpanded           = 0;                     
  my $total3PrimeContracted         = 0;                     
  my $totalExpanded                 = 0;                     
  my $totalContracted               = 0;                     
  my $totalRetrimmed                = 0;                     
  my $totalReclipped                = 0;                     
  map
  {
    my $rid = $_;
    my $disposition = ${$rh_Emit}{$rid};

    if ($disposition >= $CONSIDERED)
    {
      $aggregateClearRangeConsidered += (${$rh_OldClr}{$rid}->[1] - ${$rh_OldClr}{$rid}->[0] + 1);  # 1-based inclusive (base-based)
      $aggregateVectorFreeConsidered += (${$rh_OldClv}{$rid}->[1] - ${$rh_OldClv}{$rid}->[0] + 1);
    }

    if ($disposition >= $ASSEMBLED)
    {
      $aggregateClearRangeAssembled += (${$rh_NewClr}{$rid}->[1] - ${$rh_NewClr}{$rid}->[0] + 1);  # 1-based inclusive (base-based)
      $aggregateVectorFreeAssembled += (${$rh_NewClv}{$rid}->[1] - ${$rh_NewClv}{$rid}->[0] + 1);
    }

    if ($disposition >= $RETRIMMED)
    {
      my $expanded5prime             = (${$rh_NewClr}{$rid}->[0] < ${$rh_OldClr}{$rid}->[0])? 1 : 0; 
      my $expanded3prime             = (${$rh_NewClr}{$rid}->[1] > ${$rh_OldClr}{$rid}->[1])? 1 : 0; 
      my $contracted5prime           = (${$rh_NewClr}{$rid}->[0] > ${$rh_OldClr}{$rid}->[0])? 1 : 0; 
      my $contracted3prime           = (${$rh_NewClr}{$rid}->[1] < ${$rh_OldClr}{$rid}->[1])? 1 : 0; 
      my $expansion5prime            = ($expanded5prime)? ${$rh_OldClr}{$rid}->[0] - ${$rh_NewClr}{$rid}->[0] : 0;
      my $expansion3prime            = ($expanded3prime)? ${$rh_NewClr}{$rid}->[1] - ${$rh_OldClr}{$rid}->[1] : 0;
      my $contraction5prime          = ($contracted5prime)? ${$rh_OldClr}{$rid}->[0] - ${$rh_NewClr}{$rid}->[0] : 0;
      my $contraction3prime          = ($contracted3prime)? ${$rh_NewClr}{$rid}->[1] - ${$rh_OldClr}{$rid}->[1] : 0;
      my $expanded                   = ($expanded5prime || $expanded3prime)? 1 : 0;
      my $contracted                 = ($contracted5prime || $contracted3prime)? 1 : 0;
      my $changed                    = ($expanded || $contracted)? 1 : 0;
      my $newClearRange              = (${$rh_NewClr}{$rid}->[1] - ${$rh_NewClr}{$rid}->[0] + 1);
      my $oldClearRange              = (${$rh_OldClr}{$rid}->[1] - ${$rh_OldClr}{$rid}->[0] + 1);
      my $netChange                  = $newClearRange - $oldClearRange;    # Could be negative 

      $aggregateChange              += $netChange;
      $aggregateExpansion           += ($netChange > 0)? $netChange : 0;
      $aggregateContraction         += ($netChange < 0)? $netChange : 0;
      $aggregate5PrimeExpansion     += $expansion5prime; 
      $aggregate3PrimeExpansion     += $expansion3prime; 
      $aggregate5PrimeContraction   += $contraction5prime; 
      $aggregate3PrimeContraction   += $contraction3prime; 

      $total5PrimeExpanded += $expanded5prime;
      $total5PrimeContracted += $contracted5prime;
      $total3PrimeExpanded += $expanded3prime;
      $total3PrimeContracted += $contracted3prime;
      $totalExpanded += $expanded;
      $totalContracted += $contracted;
      $totalRetrimmed++;
    }

    if ($disposition >= $RECLIPPED)
    {
      $totalReclipped++;
    }
 
  } (keys %{$rh_Emit});
 
  # Compute means 
  my $meanClearRangeConsidered = ($nconsidered > 0)? $aggregateClearRangeConsidered / $totalReadsConsidered : 0;
  my $meanVectorFreeConsidered = ($nconsidered > 0)? $aggregateVectorFreeConsidered / $totalReadsConsidered : 0;
  my $meanClearRangeAssembled = ($nassembled > 0)? $aggregateClearRangeAssembled / $totalReadsAssembled : 0;
  my $meanVectorFreeAssembled = ($nassembled > 0)? $aggregateVectorFreeAssembled / $totalReadsAssembled : 0;
  my $meanChange                     = ($nchanged > 0)? $aggregateChange / $totalRetrimmed: 0;
  my $meanAggregateExpansion         = ($nchanged > 0)? $aggregateExpansion / $totalRetrimmed: 0;
  my $meanAggregateContraction       = ($nchanged > 0)? $aggregateContraction / $totalRetrimmed: 0;
  my $meanAggregate5PrimeExpansion   = ($nchanged > 0)? $aggregate5PrimeExpansion / $totalRetrimmed: 0;
  my $meanAggregate5PrimeContraction = ($nchanged > 0)? $aggregate5PrimeContraction / $totalRetrimmed: 0;
  my $meanAggregate3PrimeExpansion   = ($nchanged > 0)? $aggregate3PrimeExpansion / $totalRetrimmed: 0;
  my $meanAggregate3PrimeContraction = ($nchanged > 0)? $aggregate3PrimeContraction / $totalRetrimmed: 0;
    
  # Emit outputs
 
  my $of = new IO::File("> $outfilename") or $tf->bail("Cannot open .metrics file \'$outfilename\' ($!)"); 

  $of->print("readset_$lid\n");
  $of->print("\ttotalReadsConsidered = $totalReadsConsidered\n");
  $of->print("\ttotalReadsAssembled = $totalReadsAssembled\n");
  $of->print("\ttotalReadsChanged = $totalReadsChanged\n");
  $of->print("\ttotal5PrimeExpanded = $total5PrimeExpanded\n");
  $of->print("\ttotal5PrimeContracted = $total5PrimeContracted\n");
  $of->print("\ttotal3PrimeExpanded = $total3PrimeExpanded\n");
  $of->print("\ttotal3PrimeContracted = $total3PrimeContracted\n");
  $of->print("\ttotalExpanded = $totalExpanded\n");
  $of->print("\ttotalContracted = $totalContracted\n");
  $of->print("\ttotalRetrimmed = $totalRetrimmed\n");
  $of->print("\ttotalReclipped = $totalReclipped\n");
  $of->print("\tmeanClearRangeConsidered = $meanClearRangeConsidered\n");
  $of->print("\tmeanVectorFreeConsidered = $meanVectorFreeConsidered\n");
  $of->print("\tmeanClearRangeAssembled = $meanClearRangeAssembled\n");
  $of->print("\tmeanVectorFreeAssembled = $meanVectorFreeAssembled\n");
  $of->print("\tmeanChange             = $meanChange\n");
  $of->print("\tmeanAggregateExpansion = $meanAggregateExpansion\n");
  $of->print("\tmeanAggregateContraction = $meanAggregateContraction\n");
  $of->print("\tmeanAggregate5PrimeExpansion = $meanAggregate5PrimeExpansion\n");
  $of->print("\tmeanAggregate5PrimeContraction = $meanAggregate5PrimeContraction\n");
  $of->print("\tmeanAggregate3PrimeExpansion = $meanAggregate3PrimeExpansion\n");
  $of->print("\tmeanAggregate3PrimeContraction = $meanAggregate3PrimeContraction\n");
  
  $of->close() or $tf->bail("Failed to close output .metrics file \'$outfilename\' ($!)");
}


sub toList($$$$$$)
{
  my ($outprefix, $rh_Emit, $rh_OldClv, $rh_OldClr, $rh_NewClv, $rh_NewClr) = @_;
  my $outfilename;
  my @ClrContracted = ();
  my @ClrExpanded = ();
  my @ClrExchanged = ();
  my @ClrUnchanged = ();

  $outfilename = $outprefix . ".reclip.clr.contracted";
  printList($outfilename, \@ClrContracted);
  $outfilename = $outprefix . ".reclip.clr.expanded";
  printList($outfilename, \@ClrExpanded);
  $outfilename = $outprefix . ".reclip.clr.exchanged";
  printList($outfilename, \@ClrExchanged);
  $outfilename = $outprefix . ".reclip.clr.unchanged";
  printList($outfilename, \@ClrUnchanged);
}

sub getFilterSet($$)
{
  my ($infile, $rh_FilterSet) = @_;

  my $if = new IO::File("< $infile") or $tf->bail("Cannot open filter file \'$infile\' ($!)");

  my $nfiltered = 0;
  while (my $line = $if->getline())
  {
    chomp $line;
    my $rid = $line; 
    ${$rh_FilterSet}{$rid} = 1;
    ++$nfiltered;
  }

  $if->close();

  return $nfiltered;
}

# ============================================== MAIN =============================================
#
MAIN:
{    
  my %options = ();
  $options{list}       = undef;   # Whether to list out clip tables
  $options{metrics}    = undef;   # Whether to print metrics file
  $options{outprefix}  = undef;   # The prefix for output file names
  $options{accfile}    = undef;   # The name of the input acc2src lookup
  $options{featfile}   = undef;   # The name of the input .seq.features file 
  $options{resultfile} = undef;   # The name of the assembly result file
  $options{tasmfile}   = undef;   # The name of the input .tasm file
  $options{asmfile}    = undef;   # The name of the input .asm file
  $options{clrfile}    = undef;   # The name of the input .csv file
  $options{user_readref} = undef;       # User supplied read ref mode 
  $options{readref}      = undef;       # read ref mode 
  $options{lid}        = $ALL;    # Which library of sequencing to process
  $options{filterfile} = undef;   # Which subset of reads to report on, if any

  # These tables need to be populated regardless of input filtering modes.
  # 
  my %FilterSet = (); # List of subset of read names/ids to consider.
  my %Emit = ();      # List of read records and their processing status 
  my %Clb = ();       # Lookup of existing CLB range for reads as specified by input .seq.features file
  my %NewClb = ();    # Lookup of new CLB range for reads as specified by input .clb file
  my %Clr = ();       # Lookup of existing CLR range for reads from input .seq.features file
  my %NewClr = ();    # Lookup of new CLR range for reads as specified by input assembly file or .clr file
  my %Clv = ();       # Lookup of existing CLV range for reads from input .seq.features file
  my %NewClv = ();    # Lookup of new CLV range for reads as recomputed by reclip or .clv file 
  my %Rid2Alias = (); # Lookup of REad id to read alias (seq_name)
  my %Alias2Rid = (); # Lookup of read alias (seq_name) to read ID (euid)
  my @ClrContracted = ();  # Reads whose CLV remained unchanged but whose CLR contracted in assembly
  my @ClrExpanded  = ();   # Reads whose CLV remained unchanged but whose CLR expanded  in assembly
  my @ClvExchanged = ();   # Reads whose vector free interval was reset to a non-intersecting interval
  my @ClvUnchanged = ();   # Reads whose vector free interval was reset to a non-intersecting interval
  
  # Configure TIGR Foundation
  $tf->setHelpInfo($HELPTEXT.$MOREHELP);
  $tf->setUsageInfo($HELPTEXT);
  $tf->setVersionInfo($VERSION);
  $tf->addDependInfo(@DEPENDS);
  
  # validate input parameters
  my $output_options = undef;
  my $result = $tf->TIGR_GetOptions
               (
                "list!"    =>  \$options{list},
                "metrics!" =>  \$options{metrics},
                'o:s'      =>  \$options{outprefix},
                'l:s'      =>  \$options{lid},
                'ref:s'    =>  \$options{user_readref},
                'N=s'      =>  \$options{filterfile},
               );
  $tf->printUsageInfoAndExit() if (!$result || ! defined $ARGV[0]);

  # Update debug level if it's specified by the user
  my $debug_supplied = $tf->getDebugLevel();
  if (defined $debug_supplied)
  {
    $tf->setDebugLevel($debug_supplied);
    $debug = $debug_supplied;
  }
  else
  {
    $tf->setDebugLevel($debug);
  }
  setFoundation($tf);   # this is a AsmLib method that should really
                        # be part of its constructor
 
  $tf->bail("Read reference modes are : " . join (",", (keys %REF_MODE)) )
    if (defined $options{user_readref}  &&  (! exists $REF_MODE{$options{user_readref}}) );
 
  # Obtain input files       
  my @infiles = ();
  for (my $i=0; $i <= $#ARGV; $i++)
  {
    push @infiles, $ARGV[$i]; 
    $tf->bail("Cannot access input file \'$ARGV[$i]\' ($!)") if (! -r $ARGV[$i]);
    my ($prefix, $path, $suffix) = fileparse($ARGV[$i], @INPUT_SUFFIXES);
    $options{featfile} = "$prefix$suffix" if ($suffix eq "seq.features");
    $options{asmfile}  = "$prefix$suffix" if ($suffix eq "asm");
    $options{tasmfile} = "$prefix$suffix" if ($suffix eq "tasm");
    $options{clrfile}  = "$prefix$suffix" if ($suffix eq "clr");
    if (defined $options{asmfile})
    {
      $options{resultfile} = $options{asmfile};
      $options{readref} = (defined $options{user_readref})? $REF_MODE{$options{user_readref}} : $REF_BY_RID; 
    }
    elsif (defined $options{tasmfile})
    {
      $options{resultfile} = $options{tasmfile};
      $options{readref} = (defined $options{user_readref})? $REF_MODE{$options{user_readref}} : $REF_BY_ALIAS; 
    }
    elsif (defined $options{clrfile})
    {
      $options{resultfile} = $options{clrfile};
      $options{readref} = (defined $options{user_readref})? $REF_MODE{$options{user_readref}} : $REF_BY_ALIAS; 
    }
    else
    {
      $options{resultfile} = undef; 
    }
  }

  # Set .warnings filename
  if (defined $options{outprefix} && $options{outprefix} ne "")
  {
    unlink $warningsfile;
    $warningsfile = $options{outprefix} . "warnings";
    setWarnFile($warningsfile);
  }

  # Verify existence of needed files
  $tf->bail("Need a prefix.seq.features file produced by pullfrag.") if (! defined $options{featfile});  
  $tf->bail("Cannot access input file \'$options{featfile}\' ($!).") if (! -r $options{featfile});  
  $tf->bail("Need an assembly result file.") if (! defined $options{resultfile});  
  $tf->bail("Cannot access input file \'$options{resultfile}\' ($!).") if (! -r $options{resultfile});  

  # Get the filter set if there is one, undef otherwise
  my $nfiltered = 0;
  if ($options{filterfile})
  {
    progress("Using subset of reads specified in filter file \'$options{filterfile}\'");
    $tf->bail("Cannot access filter file \'$options{filterfile}\' ($!).") if (! -r $options{filterfile});  
    $nfiltered = getFilterSet($options{filterfile}, \%FilterSet);
  }

  # Phase 1 - Build sequence features from input sequence data files 
  #
  progress("Phase 1: Reading input sequencing files...");
  my $nconsidered = getClipsFromSeqFeatures($options{featfile}, \%Clv, \%Clb, \%Clr, 
                                            $options{filterfile}, \%FilterSet, \%Rid2Alias, \%Alias2Rid, \%Emit); 
  $tf->bail("Not considering any reads, terminating.") if ($nconsidered == 0);
  progress("$nconsidered reads from library \'$options{lid}\' are being considered for reclipping.");
  
  # Phase 2 - Scan input assembly or retrim files to obtain new clip points
  #           This produces the "assembled set" of reads.  This can be masked off
  #           with the -N option to focus the reads on a particular subset.
  #
  progress("Phase 2: Reading input assembly files..."); 
  my $nassembled = 0;
  if ($options{asmfile})
  {
    $nassembled = getClipsFromAsm($options{asmfile}, \%NewClr, \%Emit, \%Rid2Alias);
  } 
  elsif ($options{tasmfile})
  {
    $nassembled = getClipsFromTasm($options{tasmfile}, \%NewClr, \%Emit, \%Alias2Rid);
  }
  elsif ($options{clrfile})
  {
    # A .clr file is required but the other two are optional
    $nassembled = getClipsFromCsv($options{clrfile}, \%NewClr, \%Emit, \%Rid2Alias);
    getClipsFromCsv($options{clvfile}, \%NewClv, \%Emit, \%Rid2Alias) if (defined $options{clvfile});
    getClipsFromCsv($options{clbfile}, \%NewClb, \%Emit, \%Rid2Alias) if (defined $options{clbfile});
  }
  else
  {
    $tf->bail("No input clip files selected...");
  }
  map
  {
    $Emit{$_} = $ASSEMBLED;  # THese reads survived assembly
  } keys %NewClr;
  progress("$nassembled reads of the considered set were assembled.");
 
  # Phase 3 - Determine which reads actually need to be changed as a result of
  #           the assembly or reprocessing.  This produces the "clipped set" of reads.
  #
  progress("Phase 3: reclipping ..."); 
  my $nchanged = 0; 
  foreach my $rid (keys %Emit)
  {
    next if ($Emit{$rid} < $ASSEMBLED);    # No use checking unassembled reads 
    my $old_clrl = $Clr{$rid}->[0];
    my $old_clrr = $Clr{$rid}->[1];
    my $new_clrl = $NewClr{$rid}->[0];
    my $new_clrr = $NewClr{$rid}->[1];
    if ($new_clrl != $old_clrl  || $new_clrr != $old_clrr)
    {
      my $old_clvl = $Clv{$rid}->[0];
      my $old_clvr = $Clv{$rid}->[1];
      my ($new_clvl, $new_clvr, $disposition) = computeNewCLV($new_clrl, $new_clrr, $old_clvl, $old_clvr);
      $NewClv{$rid}->[0] = $new_clvl; 
      $NewClv{$rid}->[1] = $new_clvr;
      $NewClb{$rid}->[0] = $Clb{$rid}->[0];   # copy unchanged
      $NewClb{$rid}->[1] = $Clb{$rid}->[1];   # copy unchanged
      $Emit{$rid} = $disposition; 
      $nchanged++;
    }
  }
  progress("$nchanged reads were changed.");

  progress("Phase 4: Emitting output files...");                 
  if (defined $options{outprefix})
  {
    my $outfilename = $options{outprefix} .  ".bcp";
    toBCP($outfilename, \%Emit, \%NewClb, \%NewClv, \%NewClr, \%Rid2Alias);
  }
  if (defined $options{list})
  {
    toList($options{outprefix}, \%Emit, \%Clv, \%Clr, \%NewClv, \%NewClr);
  }
  if (defined $options{metrics})
  {
    my $outfilename = $options{outprefix} . ".metrics";
    progress("Writing .metrics file \'$outfilename\' ...");
    toMetrics($outfilename, \%Emit, \%Clv, \%Clr, \%NewClv, \%NewClr, $options{lid}, $nconsidered, $nassembled, $nchanged);
  }
 
  progress("Done.");
  exit 0;
}

