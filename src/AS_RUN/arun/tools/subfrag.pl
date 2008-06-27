#!/usr/local/bin/perl -w
# $Id: subfrag.pl,v 1.2 2008-06-27 06:29:19 brianwalenz Exp $
# subfrag - Extract a filter of reads (and possibly their mates) from a
#           set of .frg files.
#
# Written by Martin Shumway based on an earlier script by Daniela Puiu.
#

use strict;
use warnings;
use File::Basename;

# TIGR Modules
use TIGR::Foundation;
use TIGR::AsmLib;

my $tf = new TIGR::Foundation;
my $PRG = $tf->getProgramInfo('name');
my $VERSION="1.54";
my @DEPENDS=("TIGR::Foundation","TIGR::AsmLib");
my $SORT = "/usr/bin/sort";
push @DEPENDS, $SORT;

my $HELPTEXT = qq~
Extract fragments from a CA .frg file.

subfrag  <frg1 frg2 ...>  [options]

  frgs   One or more .frg files as constructed by pullfrag and related tools.
         File list can be resolved by the shell, for example my*.frg

  options:
    -N <file>    Optional file with the read names to be extracted (one/line).
    -s read[,read]...  Comma separated list (no spaces) of read ids.
    --mates      Extract mates as well
    --[no]links  Report mate links (default).  --nolinks suppresses links.
    -3 <K>       Trim back K bases from the 3' end of the fragment sequence.
                 The K value can be negative as well (to extend the read).
    -5 <K>       Trim back K bases from the 5' end of the fragment sequence.
                 The K value can be negative as well (to extend the read).
    -1 <file>    Retrim specified reads from retrim file (1-based inclusive)
    -0 <file>    Retrim specified reads from retrim file (0-based exclusive)
    -o <files>   Write output to file instead of the console.  Multiple
                 output option files can be specified (comma separated):
                 <prefix>.frg : Single output .frg file
                 <prefix>.seq : sequence of entire read (TIGR form)
                 <prefix>.seqs : List of reads by name (or accession if no name)
                 <prefix>.qual: qualities of entire read (phred style)
                 <prefix>.clr : Clear range (rid clrl clrr) CA form.
                 <prefix>.clr.fasta : sequence of clear range of read
                 <prefix>.acc : Lookup between acc,src fields, as available
                 <prefix>.mates : Mate pairs encountered
                 <prefix>.bambus.mates : Mate pairs encountered (bambus format)
                 <prefix>.ftab : rid/clr.fasta in single tab delimited line

Given one or more input .frg files and an optional filter file, subfrag extracts
a subset of fragments and builds a new .frg file with associated DST, LKG and
FRG records.

If no filter option is specified then all records are emitted on output.
Specifying multiple files on the command line has the effect of merging the
FRG records prior to filtering.  Duplicate records are omitted.

The read set can be decimated using the -N or -s options.  If no qualifier
is supplied, then entire input read set is used (default).  Input reads can be
specified by accession or read name (if assigned).  When filtering you can
specify that mates be included in the output set (--mates).

Some trimming operations can be performed on the output set.  3' end 5' end trim
can be performed on all reads using the -3 option, respectively.  Also, a mask
file with new trim coordinates can be applied to specified reads either in
classical form or 0-based CA form.  The newly trimmed sequence must meet the CA
criteria for fragment length.  Trims that don't are not performed and a message
is written to the .warnings file.

subfrag is useful for performing in-situ operations on .frg files without going
through a database intermediary.  Additional output conversions can be specified
using the -o option.  Otherwise, a new .frg file is generated on stdout.
~;

my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;

# =============================== Constants ================================

# Input/Output parameters
my $NAME = "src";
my $EUID = "acc";
my %SOURCE_MODES = ("name" => $NAME, "euid" => $EUID);
my @OUTPUT_SUFFIXES =
(
  "frg",
  "seq",
  "clr.fasta",
  "clr",
  "qual",
  "seqs",
  "acc",
  "bambus.mates",
  "mates",
  "ftab",
);
my %SUPPORTED_OUTPUTS = ();
map { $SUPPORTED_OUTPUTS{$_} = 1; } @OUTPUT_SUFFIXES;
my $QUAL_LEN    = 17;      # number of quality values in one line of the
                           # .qual file

# Operational parameters
my $debug = 0;
my $LOW = 1;
my $MEDIUM = 2;
my $HIGH = 3;
my $PROGRESS_FACTOR = 100000;
my $outprefix = "tmp";
my $warningsfile = "$outprefix.warnings";

# Trimming parameters
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

# getSid - Return the normalized sequence id (sid) given the frg record of
# acc (sequence accession) and optional src (sequence name).
#
sub getSid($)
{
  my ($rh_fields) = @_;
  my $sid = undef;

  my %fields = %$rh_fields;
  my $acc  = $fields{acc};
  if (defined $fields{src}  &&  $fields{src} ne "")
  {
    $sid  = $fields{src};
    chop $sid;
  }
  else
  {
    $sid = $acc;
  }

  return $sid;
}

# setClr - Procedure to reset clip points
#   Rule 1:  Apply 3' trimming to all reads, if specified
#   Rule 2:  Apply 5' trimming to all reads, if specified
#   Rule 3:  Apply custom trimming to specified reads, overriding
#          any other trimming.
#   The Trim for a read is canceled if during any sequence of trimming
#   operations one of the following constraints is violated:
#     a. CLR is too short for assembly in CA
#     b. CLR is too long for assembly in CA
#     c. CLR intersects boundary of raw sequence
#   In this case the frag is emitted anyway with the unmodified clr points
#   Trim warnings are delivered in the outprefix.warnings file
#

sub setClr($$$$$$$)
{
  my ($rh_fields, $trim3Prime, $trim5Prime, $trimca, $trimta, $rh_Clr, $rh_NewClr) = @_;
  my %fields = %{$rh_fields};

  my $sid = getSid(\%fields);

  my ($clrl, $clrr) = split /,/,$fields{clr};     # Initialize
  my $clrl_orig = $clrl;
  my $clrr_orig = $clrr;
  my $sequence = $fields{seq};
  $sequence =~ s/\s//g;
  my $seqlen = length($sequence);

  if ($trim5Prime)
  {
    $clrl = $clrl + $trim5Prime;
  }
  if ($trim3Prime)
  {
    $clrr = $clrr - $trim3Prime;
  }

  if ($trimca || $trimta)
  {
    my $sid = getSid(\%fields);
    if (exists ${$rh_NewClr}{$sid})
    {
      $clrl = ${$rh_NewClr}{$sid}->[0];
      $clrr = ${$rh_NewClr}{$sid}->[1];
    }
  }

  # Now check new settings for validity.
  my $invalid = 0;
  $invalid = $FRG_TOO_SHORT if ($clrr  - $clrl  < $MINLENGTH);
  $invalid = $FRG_TOO_LONG  if ($clrr  - $clrl  > $MAXLENGTH);
  $invalid = $FRG_OUT_RANGE if ($clrl  > $seqlen-1  ||  $clrr  > $seqlen);

  if ($invalid)
  {
    ${$rh_Clr}{$sid} = \[$clrl_orig, $clrr_orig];
    logwarning("Fragment $sid has $FRG_VALID_MESSAGE{$invalid}.  " .
               "Ignoring attempt to adjust clear range to [$clrl, $clrr).", 1);
  }
  else
  {
    ${$rh_Clr}{$sid} = \[$clrl, $clrr];
  }
}

# ============================================== MAIN =============================================
#
MAIN:
{
  my %options = ();
  $options{trim3Prime} = undef;
  $options{trim5Prime} = undef;
  $options{trimca}     = undef;
  $options{trimta}     = undef;
  $options{filter} = 0;
  $options{filterfile} = undef;
  $options{filterlist} = undef;
  $options{links}      = 1;

  # These tables need to be populated regardless of input filtering modes.
  #
  my %Mates = ();     # List of all mate pairs discovered in input
  my %Emit  = ();     # Track which FRG records to be emitted by acc
  my %Src = ();       # List of input FRG records indexed by acc, with optional src
  my %FilterSet = (); # List of reads (by acc or by src) to emit as a subset
  my %Dst = ();       # List of input DST records indexed by acc
  my %Lkg = ();       # List of input LKG records indexed by acc member
  my %Lkg2Dst = ();   # Lookup table of Dst acc given Lkg record
  my %Clr = ();       # Lookup of existin clr range for reads from existing .frg records
  my %NewClr = ();    # Lookup of new clr range for reads as specified by retrim file

  # Configure TIGR Foundation
  $tf->setHelpInfo($HELPTEXT.$MOREHELP);
  $tf->setUsageInfo($HELPTEXT);
  $tf->setVersionInfo($VERSION);
  $tf->addDependInfo(@DEPENDS);

  # validate input parameters
  my $output_options = undef;
  my $result = $tf->TIGR_GetOptions
               (
                'N:s'     =>  \$options{filterfile},
                's:s'     =>  \$options{filterlist},
                "mates!"  =>  \$options{mates},
                "links!"  =>  \$options{links},
                '3=i'     =>  \$options{trim3Prime},
                '5=i'     =>  \$options{trim5Prime},
                '0=s'     =>  \$options{trimca},
                '1=s'     =>  \$options{trimta},
                'o:s'     =>  \$output_options,
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

  my @infiles = ();
  for (my $i=0; $i <= $#ARGV; $i++)
  {
    push @infiles, $ARGV[$i];
    $tf->bail("Cannot access input file \'$ARGV[$i]\' ($!)") if (! -r $ARGV[$i]);
    my ($newoutprefix, $path, $suffix) = fileparse($ARGV[$i], @OUTPUT_SUFFIXES);
    if (defined $newoutprefix && $newoutprefix ne "")
    {
      $warningsfile = $newoutprefix . "warnings";
    }
    $options{warningsfile} = $warningsfile;
    $options{outprefix} = $newoutprefix;
  }
  unlink $warningsfile;
  setWarnFile($warningsfile);

  # Check whether to include LKG output
  progress("Mate pair information suppressed by request.") if (! $options{links});

  # Get optional list of filter of sequences to output.  If this
  # option is not supplied then the whole input set is used.
  #
  if (defined $options{filterfile})
  {
    $options{filter} = 1;

    #get read names from input file
    my $ff = new IO::File("< $options{filterfile}")
      or $tf->bail("Failed to open input filter file \'$options{filterfile}\' ($!)");
    while (my $line = $ff->getline())
    {
      chop $line;
      $line =~ s/\s+//g;
      next if ($line =~/^#/);
      $FilterSet{$line} = 1;
    }
    $ff->close();
  }
  elsif (defined $options{filterlist})
  {
    $options{filter} = 1;
    my @f = split ",",$options{filterlist};
    foreach my $read (@f)
    {
      $FilterSet{$read} = 1;
    }
  }
  else
  {
    progress("No filter specified, using all records.");
  }

  # Get possible trim spec from file.  Build table.
  if (defined $options{trimca})
  {
    my $trimfilename = "$options{trimca}";
    my $trimf = new IO::File("< $trimfilename") or $tf->bail("Could not open retrim file \'$trimfilename\' ($!)");
    while (my $line = $trimf->getline())
    {
      chop $line;
      my @f = split /\s+/,$line;
      my $sid = $f[0];
      my $clrl = $f[1];
      my $clrr = $f[2];
      $NewClr{$sid} = [$clrl,$clrr];
    }
    $trimf->close();
  }
  elsif (defined $options{trimta})
  {
    my $trimfilename = "$options{trimta}";
    my $trimf = new IO::File("< $trimfilename") or $tf->bail("Could not open retrim file \'$trimfilename\' ($!)");
    while (my $line = $trimf->getline())
    {
      chop $line;
      my @f = split /\s+/,$line;
      my $sid = $f[0];
      my $clrl = $f[1] - 1;  # normalize to CA format
      my $clrr = $f[2];
      $NewClr{$sid} = [$clrl,$clrr];
    }
    $trimf->close();
  }

  # Get output options
  my %Output_options = ();
  my $hasNamedOutputStream = 0;
  if (defined $output_options)
  {
    my @outputs = split /,/,$output_options;
    map
    {
      my ($name, $path, $suffix) = fileparse($_, @OUTPUT_SUFFIXES);
      $tf->bail("Unsupported output option in -o specification: \'$_\'") if (! exists $SUPPORTED_OUTPUTS{$suffix});
      $Output_options{$suffix} = ($path eq "./")? "$name$suffix" : "$path/$name$suffix";
      $hasNamedOutputStream = 1;
    } @outputs;
    $options{frgout} = (exists $Output_options{frg})? $Output_options{frg} : undef;
    $options{seqout} = (exists $Output_options{seq})? $Output_options{seq} : undef;
    $options{seqsout} = (exists $Output_options{seqs})? $Output_options{seqs} : undef;
    $options{qualout} = (exists $Output_options{qual})? $Output_options{qual} : undef;
    $options{clrout} = (exists $Output_options{clr})? $Output_options{clr} : undef;
    $options{fastaout} = (exists $Output_options{"clr.fasta"})? $Output_options{"clr.fasta"} : undef;
    $options{accout} = (exists $Output_options{acc})? $Output_options{acc} : undef;
    $options{matesout} = (exists $Output_options{mates})? $Output_options{mates} : undef;
    $options{bambusmatesout} = (exists $Output_options{"bambus.mates"})? $Output_options{"bambus.mates"} : undef;
    $options{ftabout} = (exists $Output_options{ftab})? $Output_options{ftab} : undef;
  }

  printTable($options{outprefix} . "options", \%options) if ($debug > 0);

  # PASS 1: Build lookup tables of objects seen.
  #
  # Scan the .frg files and build various lookup tables.
  # This has to be done using the accession (acc) values.
  # This data needed regardless of filtering options.
  #
  progress("Phase 1 : Scan input .frg files ...");

  my $trim3Prime = (exists $options{trim3Prime})? $options{trim3Prime} : undef;
  my $trim5Prime = (exists $options{trim5Prime})? $options{trim5Prime} : undef;
  my $trimca     = (exists $options{trimca})?     $options{trimca}     : undef;
  my $trimta     = (exists $options{trimta})?     $options{trimta}     : undef;

  foreach my $infile (@infiles)
  {
    my $if = new IO::File("< $infile") or $tf->bail("Cannot open input .frg file \'$infile\' ($!)");
    progress("Phase 1 : Scanning input file \'$infile\'...");

    my $nrecords = 0;
    while (my $rec = getCARecord($if))
    {
        my ($type, $rh_fields, $recs) = parseCARecord($rec);
        my %fields = %{$rh_fields};
        if ($type eq "FRG")
        {
          my $acc = $fields{acc};
          my $sid = getSid(\%fields);
          $Src{$acc} = $sid;
          setClr(\%fields, $trim3Prime, $trim5Prime, $trimca, $trimta, \%Clr, \%NewClr);
        }
        elsif ($type eq "LKG")
        {
          my $fg1 = $fields{fg1};
          my $fg2 = $fields{fg2};
          my $dst = $fields{dst};

          $Mates{$fg1} = $fg2;
          $Mates{$fg2} = $fg1;
          $Lkg{$fg1} = 1;       # track link record based on at least one constituent accession
          $Lkg2Dst{$fg1} = $dst;  # track canonical mate pair as the one listed first
        }

      ++$nrecords;
      progress("Phase 1: $nrecords records scanned.") if ($nrecords > 1 && $nrecords % $PROGRESS_FACTOR == 1);
    }

    $if->close();
  }

  # Setup output streams
  my $of = undef;
  my $sf = undef;
  my $nf = undef;
  my $cf = undef;
  my $qf = undef;
  my $ff = undef;
  my $mf = undef;
  my $bf = undef;
  my $af = undef;
  if ($hasNamedOutputStream)
  {
    if (defined $options{frgout})
    {
      my $outfilename = $options{frgout};
      $of = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{seqout})
    {
      my $outfilename = $options{seqout};
      $sf = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{seqsout})
    {
      my $outfilename = $options{seqsout};
      $nf = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{qualout})
    {
      my $outfilename = $options{qualout};
      $qf = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{clrout})
    {
      my $outfilename = $options{clrout};
      $cf = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{fastaout})
    {
      my $outfilename = $options{fastaout};
      $ff = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{matesout})
    {
      my $outfilename = $options{matesout};
      $mf = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{bambusmatesout})
    {
      my $outfilename = $options{bambusmatesout};
      $bf = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
    if (defined $options{ftabout})
    {
      my $outfilename = $options{ftabout};
      $af = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    }
  }
  else
  {
    $of = new IO::Handle;
    $of->fdopen(fileno(STDOUT), "w");
  }


  # PASS 2: Emit FRG records
  #
  # Stream out the input file such that each FRG record is read and emitted
  # with possible modifications depending on execution mode.
  #
  progress("Phase 2 : Emit FRG records...");
  my $firstfile = 1;
  foreach my $infile (@infiles)
  {
    my $if = new IO::File("< $infile") or $tf->bail("Cannot open input .frg file \'$infile\' ($!)");
    progress("Phase 2 : Scanning Scanning input file \'$infile\'...");

    # Stream through the input .frg files a second time and emit
    # the dataset depending on filtering options, and optionally
    # perform terminal processing on the records.
    #
    my $nrecords = 0;
    while (my $rec = getCARecord($if))
    {
      my $isPrint = 0;
      my($id, $rh_fields, $recs) = parseCARecord($rec);
      my %fields=%{$rh_fields};

      if ($id eq "FRG")
      {
        if ( $options{filter} )
        {
          my $acc  = $fields{acc};
          my $sid = getSid(\%fields);
          # Check to see first whether the record is in the filter set,
          # then whether it is a mate of a record in the filter set and
          # we are to draw the mates.
          if (exists $FilterSet{$sid})
          {
            if (! exists $Emit{$acc})
            {
              $isPrint = 1;
              $Emit{$acc} = $acc;        # track the set of reads to emit by acc
            }
            else
            {
              $isPrint = 0;     # suppress duplicate records
            }
          }
          elsif ($options{mates})
          {
            if (exists $Mates{$acc})
            {
              my $mate_acc = $Mates{$acc};
              my $mate_src = $Src{$mate_acc};
              if ( (exists $FilterSet{$mate_acc})  ||  (exists $FilterSet{$mate_src}) )
              {
                if (! exists $Emit{$acc})
                {
                  $isPrint = 1;
                  $Emit{$acc} = $acc;        # track the set of reads to emit by acc
                }
                else
                {
                  $isPrint = 0;     # suppress duplicate records
                }
              }
            }
          }
        }
        else
        {
          my $acc  = $fields{acc};
          if (! exists $Emit{$acc})
          {
            $isPrint = 1;
            $Emit{$acc} = $acc;        # track the set of reads to emit by acc
          }
          else
          {
            $isPrint = 0;     # suppress duplicate records
          }
        }
      }
      elsif($id eq "LKG")
      {
        if ($options{filter})
        {
          # Determine which LKG records to keep based on whether the Emit set
          # contains both elements of the LKG object.  Note that the LKG record
          # must come after the constituent FRGs, as is required in the
          # specification.  This lookup is by acc only.
          #
          my $fg1 = $fields{fg1};
          my $fg2 = $fields{fg2};
          if ( (exists $Emit{$Mates{$fg1}})  &&  (exists $Emit{$Mates{$fg2}}) )
          {
            if (exists $Lkg{$fg1} && $Lkg{$fg1} > 0)
            {
              $isPrint = 1;
              $Lkg{$fg1} = 0;   # Mark as emitted
            }
          }
        }
        else
        {
          my $fg1 = $fields{fg1};   # Use fg1 only as surrogate for entire LKG record
          if (exists $Lkg{$fg1} && $Lkg{$fg1} > 0)
          {
            $isPrint = 1;
            $Lkg{$fg1} = 0;   # Mark as emitted
          }
        }
        $isPrint = 0 if (! $options{links});   # Suppress mate pair linkages
      }
      elsif ($id eq "DST")
      {
      	my $dst = $fields{acc};
        $Dst{$dst} = (exists $Dst{$dst})? $Dst{$dst} + 1 : 1;
        # Emit only first encounter with a DST record of a given accession
        if ($Dst{$dst} == 1)
        {
          $isPrint = 1;
        }
        else
        {
          $isPrint = 0;
        }
      }
      elsif ($id eq "BAT")
      {
        $isPrint = ($firstfile)? 1 : 0;   # suppress multiple ADT records
      }
      elsif ($id eq "ADT")
      {
        if ($firstfile)
        {
          my $username = getpwuid($<);
          $rec =~ s/who:(\w+)/who:$username/o;
          my $agent = "Generated by $PRG $VERSION";
          $rec =~ s/com:\n*\n/com:\n$agent\n/o;
          $isPrint = 1;
        }
        else
        {
          $isPrint = 0;   # suppress subsequent ADT records
        }
      }
      else
      {
        $isPrint = 1;  # pass through all non-FRG and non-LKG records
      }

      if ($isPrint)
      {
        if ($id eq "FRG")
        {
          my $sid = getSid(\%fields);

          my $r = $Clr{$sid};
          my $clrl = $$r->[0];
          my $clrr = $$r->[1];
          # This hack exists because AsmLib has no CA record constructor
          $fields{clr} = "$clrl,$clrr";
          my $replacement = "clr:$fields{clr}";
          $rec =~ s/clr:(\d+),(\d+)/$replacement/o;

          if (defined $sf)
          {
            my $sequence = $fields{seq};
            my $clrl_modified = $clrl + 1;  # modify to 1-based inclusive
            $sf->print(">$sid    0 0 0   $clrl_modified $clrr\n" . uc $sequence) or $tf->bail("Failed to write to output seq stream ($!)");
          }

          if (defined $nf)
          {
            $nf->print("$sid\n") or $tf->bail("Failed to write to output .seqs stream ($!)");
          }

          if (defined $qf)
          {
            my $qualityS = $fields{qlt};
            $qualityS =~ s/\s+//g;
            my $qualities = "";
            map
            {
              $qualities .= sprintf("%02d ", ord($_) - 48);
            } (split //,$qualityS);
            $qualities =~ s/((\d\d\s+){$QUAL_LEN})/$1\n/g;
            $qualities =~ s/\s+$//g;
            $qf->print(">$sid\n$qualities\n") or $tf->bail("Failed to write to output qual stream ($!)");
          }

          if (defined $ff || defined $af)
          {
            # convert clrr to 0-based inclusive so substr can use it
            my $clrr_modified = $clrr-1;
            my $sequence = $fields{seq};
            $sequence =~ s/\s//g;
            my $clr_sequence = substr($sequence, $clrl, $clrr_modified - $clrl+1);
            my $length = length($clr_sequence);
            if (defined $af)
            {
              $af->print("$sid\t$clr_sequence\n") or $tf->bail("Failed to write to output ftab stream ($!)");
            }
            $clr_sequence =~ s/((.){60})/$1\n/g;
            chomp $clr_sequence if ($length % 60 == 0);   # don't add \n to end of full line
            if (defined $ff)
            {
              $ff->print(">$sid\n$clr_sequence\n") or $tf->bail("Failed to write to output clr.fasta stream ($!)");
            }
          }

          if (defined $cf)
          {
            $cf->print("$sid\t$clrl\t$clrr\n") or $tf->bail("Failed to write to output clr stream ($!)");
          }
        }
        elsif ($id eq "DST")
        {
          if (defined $bf)
          {
            my $acc = $fields{acc};
            my $mean = $fields{mea};
            my $stddev = $fields{std};
            my $min = ($mean - $stddev > 0)? $mean - $stddev : 1;
            my $max = $mean + $stddev;
            my $line = "library\t$acc\t$min\t$max\n";
            $bf->print($line) or $bf->bail("Failed to write to output bambus.mates stream ($!)");
          }
        }
        elsif ($id eq "LKG")
        {
          if (defined $mf || defined $bf)
          {
            my $fg1 = $fields{fg1};
            my $name1 = $Src{$fg1};
            my $fg2 = $fields{fg2};
            my $name2 = $Src{$fg2};
            if (defined $mf)
            {
              my $line = "$name1\t$name2\n";
              $mf->print($line) or $mf->bail("Failed to write to output mates stream ($!)");
            }
            if (defined $bf)
            {
              my $libid = $Lkg2Dst{$fg1};
              my $line = "$name1\t$name2\t$libid\n";
              $bf->print($line) or $bf->bail("Failed to write to output bambus.mates stream ($!)");
            }
          }
        }

        # Emit all records to main .frg stream
        if (defined $of)
        {
          $of->print($rec) or $tf->bail("Failed to write to output frg stream ($!)");
        }

        $nrecords++;
        progress("Phase 2: $nrecords records scanned.") if ($nrecords > 1 && $nrecords % $PROGRESS_FACTOR == 1);
      }
    }

    $if->close();
    $firstfile = 0;
  }

  progress("Phase 3: Emit auxiliary files...");

  # Close output stream
  if (defined $of)
  {
    $of->close() or $tf->bail("Failed to close output frg stream ($!)");
  }
  if (defined $cf)
  {
    $cf->close() or $tf->bail("Failed to close output clr stream ($!)");
  }
  if (defined $sf)
  {
    $sf->close() or $tf->bail("Failed to close output seq stream ($!)");
  }
  if (defined $ff)
  {
    $ff->close() or $tf->bail("Failed to close output clr.fasta stream ($!)");
  }
  if (defined $nf)
  {
    $nf->close() or $tf->bail("Failed to close output seqs stream ($!)");

    # Resort the output from this step.
    my $sortinfile = $options{seqsout};
    my $sortoutfile = "$sortinfile.tmp";
    my $cmd = "$SORT -u $sortinfile > $sortoutfile";
    $tf->logLocal("Running command: \'$cmd\'", $LOW);
    my $bad = $tf->runCommand($cmd);
    $tf->bail("Failed to run command \'$cmd\' ($!)") if ($bad);
    unlink ($sortinfile);
    rename($sortoutfile, $sortinfile) or $tf->bail("Could not rename file \'$sortoutfile\' ($!)");
  }
  if (defined $qf)
  {
    $qf->close() or $tf->bail("Failed to close output qual stream ($!)");
  }
  if (defined $mf)
  {
    $mf->close() or $tf->bail("Failed to close output .mates stream ($!)");
  }
  if (defined $bf)
  {
    $bf->close() or $tf->bail("Failed to close output .bambus.mates stream ($!)");
  }
  if (defined $af)
  {
    $af->close() or $tf->bail("Failed to close output .ftab stream ($!)");
  }

  if (defined $options{accout})
  {
    my $accfile = "$options{accout}";
    printTable($accfile, \%Src);
  }

  progress("Done.");
  exit 0;
}

