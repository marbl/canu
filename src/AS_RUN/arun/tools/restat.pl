#!/usr/local/bin/perl -w
# $Id: restat.pl,v 1.2 2008-06-27 06:29:19 brianwalenz Exp $
# restat - Extract a filter of reads (and possibly their mates) from a
#           set of .frg files.
#
# Written by Daniella Puiu and Martin Shumway
#

use strict;
use warnings;
use File::Basename;
use File::Copy "cp";
use FileHandle;

# TIGR Modules
use TIGR::Foundation;
use TIGR::AsmLib;

my $tf = new TIGR::Foundation;
setFoundation($tf);   # this is a AsmLib method that should really be part of its constructor
my $PRG = $tf->getProgramInfo('name');
my $VERSION="1.06";
my @DEPENDS=
(
  "TIGR::Foundation",
  "TIGR::AsmLib",
);

my $HELPTEXT = qq~
Extract library insert statistics from a CA assembly .asm file.

restat  <prefix>.asm  <insertfile>  [options]

  <prefix>.asm   Input .asm file from CA assembly output.
  <insertfile>   Path to .insert file produced by pullfrag or other program

  options:
    -o <files,[file]...>  Specify one or more output file types:
                   .clone  - old style format: cat\# mean stddev
                   .insert - rewrite the .insert file with new stats.
                             Any existing file by that name is moved aside.
                   .metrics - emit an INI file with aggregate stats for ametrics
                   .library - XML file format for genomic libraries and stats

The restat program reads MDI records from the Celera Assembler .asm file and
extracts genomic library insert size statistics so that the database can be
updated.  A .clone file, or other output files are generated to support upload.

SEE ALSO
  aloader
~;

# Maintenance notes:
# We could greatly improve on this software with the following changes:
# 1. Impleement the library data as a perl module whose constructor can read
# and existing .library file, and replace the code here with method calls.
# 2. Get rid of .clone file and change uploadInserts to read the .library file
# (this is planned as part of ECP 1016).
# 3. Get rid of the .insert file and replace it with a .library file to hold
# stats changes, and a static .extents xml file to hold insert info and library
# membership.  Then the latter file will not be needed by this program, which
# would instead read and write .library files.
# shumwaym 12/05/06

# =============================== Constants ================================

# Input/Output parameters
my @INPUT_SUFFIXES = ("asm", "insert");
my %SUPPORTED_INPUTS  = ();
map { $SUPPORTED_INPUTS{$_} = 1; } @INPUT_SUFFIXES;
my $CLONE  = "clone";
my $INSERT = "insert";
my $METRICS = "metrics";
my $LIBRARY = "library";
my @OUTPUT_SUFFIXES = ($CLONE, $INSERT, $METRICS, $LIBRARY);
my %SUPPORTED_OUTPUTS  = ();
map { $SUPPORTED_OUTPUTS{$_} = 1; } @OUTPUT_SUFFIXES;

# Operational parameters
my $debug = 0;
my $LOW = 1;
my $MEDIUM = 2;
my $HIGH = 3;
my $warningsfile = "tmp.warnings";
my $PROGRESS_FACTOR = 10000000;

# Formatting Constants
my %SECTION_ITEMS =
(
  'lid' => 0,
  'name' => 1,
  'desc' => 2,
  'stats' => 3,
  'histogram' => 4,
);
my %STAT_ITEMS =
(
  'input_mean' => 0,
  'mean' => 1,
  'input_stddev' => 2,
  'stddev' => 3,
  'max' => 4,
  'min' => 5,
  'N' => 6,
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

# Resort items in a table into a rank order as specified in the reference table.
#
sub rankOrder($$)
{
  my ($rh_query, $rh_reference) = @_;

  my @keys = keys %{$rh_query};
  my %rank2value = ();
  my $k=10000;
  map
  {
    my $tag = $_;
    my $rank = 0;
    if (exists ${$rh_reference}{$_})
    {
      $rank = ${$rh_reference}{$_};
    }
    else
    {
      $rank = $k++;
    }
    $rank2value{$rank} = $tag;
  } @keys;

  return \%rank2value;
}

# ============================================== MAIN =============================================
#
MAIN:
{
  my %options = ();
  $options{insertfile} = undef;
  $options{cloneoutfile} = undef;
  $options{insertoutfile} = undef;
  $options{libraryoutfile} = undef;
  $options{prefix} = undef;
  $options{outprefix} = undef;
  $options{removeEmpty} = 1;    # Remove from output libs with no experimental data

  # These tables need to be populated regardless of input filtering modes.
  #
  my %ResultLibrary = ();   # Tracked libraries
  my %InputLibrary = ();   # List of all mate pairs discovered in input

  # Configure TIGR Foundation
  $tf->setHelpInfo($HELPTEXT);
  $tf->setUsageInfo($HELPTEXT);
  $tf->setVersionInfo($VERSION);
  $tf->addDependInfo(@DEPENDS);

  # validate input parameters
  my $output_options = undef;
  my $result = $tf->TIGR_GetOptions
               (
                'o=s'     =>  \$output_options,
               );
  $tf->printUsageInfoAndExit() if (!$result || $#ARGV < 1);

  my ($prefix, $path, $suffix) = fileparse($ARGV[0], '.asm');
  $tf->bail("Need a .asm file") if (! defined $prefix);
  $options{asmfile} = "$path$prefix$suffix";
  $options{outprefix} = $prefix;
  $tf->bail("Cannot access input .asm file \'$options{asmfile}\' ($!)") if (! -r $options{asmfile});

  if (defined $ARGV[1])
  {
    my ($prefix, $path, $suffix) = fileparse($ARGV[1], '.insert');
    if (defined $prefix)
    {
      $options{insertfile} = ($path eq "./")? "$prefix$suffix" : "$path$prefix$suffix";
      $tf->bail("Cannot access input .insert file \'$options{insertfile}\' ($!)") if (! -r $options{insertfile});
    }
  }

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

  if (defined $options{outprefix} && $options{outprefix} ne "")
  {
    $warningsfile = "$options{outprefix}.warnings";
  }
  else
  {
    $warningsfile = "$options{prefix}.warnings";
  }
  unlink $warningsfile;
  setWarnFile($warningsfile);

  # Get output options
  my %Output_options = ();
  if (defined $output_options)
  {
    my @outputs = split /,/,$output_options;
    map
    {
      my ($name, $path, $suffix) = fileparse($_, @OUTPUT_SUFFIXES);
      $tf->bail("Unsupported output option in -o specification: \'$_\'") if (! exists $SUPPORTED_OUTPUTS{$suffix});
      $Output_options{$suffix} = ($path eq "./")? "$name$suffix" : "$path/$name$suffix";
    } @outputs;
    $options{cloneoutfile} = (exists $Output_options{clone})? $Output_options{clone} : undef;
    $options{insertoutfile} = (exists $Output_options{insert})? $Output_options{insert} : undef;
    if ($options{insertoutfile})
    {
      if ($options{insertoutfile} eq $options{insertfile})
      {
        my ($name, $path, $suffix) = fileparse($options{insertoutfile}, @OUTPUT_SUFFIXES);
        $options{insertfile} = ($path eq "./")? "$name" . "0.$suffix" : "$path/$name" . "0.$suffix";

        progress("Saving off input .insert file \'$options{insertoutfile}\' as it will be rewritten.");
        if (-l $options{insertfile})
        {
          progress("Input .insert file is a link, copying to \'$options{insertfile}\'...");
          my $fc = FileHandle->new($options{insertfile}, "r");
          cp($options{insertoutfile}, $options{insertfile})
            or $tf->bail("Cannot copy file \'$options{insertoutfile}\' to \'$options{insertfile}\' ($!)");
        }
        else
        {
          progress("Input .insert file exists, moving to \'$options{insertfile}\'...");
          my $result = rename($options{insertoutfile}, $options{insertfile});
          $tf->bail("Cannot rename file \'$options{insertoutfile}\' to \'$options{insertfile}\' ($!)") if ($result != 1);
        }
      }
    }
    $options{metricsfile} = (exists $Output_options{metrics})? $Output_options{metrics} : undef;
    $options{libraryoutfile} = (exists $Output_options{library})? $Output_options{library} : undef;
  }

  # PASS 1: Obtain updated MDI records from .asm file
  #
  progress("Scanning input .asm file \'$options{asmfile}\'...");
  {
    my $if = new IO::File("< $options{asmfile}") or $tf->bail("Cannot open input .asm file \'$options{asmfile}\' ($!)");

    my $nrecords = 0;
    while (my $rec = getCARecord($if))
    {
        my ($type, $rh_fields, $recs) = parseCARecord($rec);
        my %fields = %{$rh_fields};
        if ($type eq "MDI")
        {
          my $ref = $fields{ref};
          $ref =~ m/\((\d+),(\d+)\)/o;
          my $acc = $1;
          $ResultLibrary{$acc}{name}  = "name?";         # bind this later from the .insert file
          $ResultLibrary{$acc}{stats}{mean} = $fields{mea};
          $ResultLibrary{$acc}{stats}{stddev}= $fields{std};
          my $max =  $fields{max};
          my $min =  $fields{min};
          $ResultLibrary{$acc}{stats}{max}  = $max;
          $ResultLibrary{$acc}{stats}{min}  = $min;
          my $buckets = $fields{buc};
          $ResultLibrary{$acc}{hist}{buc}  = $buckets;
          my @histogram = split /\n/,$fields{his};
          $ResultLibrary{$acc}{hist}{values}  = \@histogram;
          my $N = 0;
          map { $N += ($_ =~ m/\d+/)? $_ : 0; } @histogram;
          $ResultLibrary{$acc}{stats}{items}    = $N;
          $ResultLibrary{$acc}{hist}{bin}  = ($buckets > 0)? ($max - $min) / $buckets : $max - $min;
        }

      ++$nrecords;
      progress("$nrecords .asm records scanned.") if ($nrecords > 1 && $nrecords % $PROGRESS_FACTOR == 1);
    }

    $if->close();
  }

  # PASS 2: Obtain existing Genomic library records from .insert file
  #
  progress("Scanning input .insert file \'$options{insertfile}\'...");
  my $if = new IO::File("< $options{insertfile}")
      or $tf->bail("Cannot open input .insert file \'$options{insertfile}\' ($!)");

  my $nrecords = 0;
  while (my $line = $if->getline())
  {
      chomp $line;
      # I wish we had a proper module to interact with this file type
      if ($line =~ m/\s+<GenomicLibrary Id="(\w+)" acc="(\d+)" mea="(\d+).(\d+)" std_dev="(\d+).(\d+)">/o)
      {
        my $name = $1;
        my $acc = $2;
        my $oldmean = sprintf("%s.%s", $3, $4) + 0.0;
        my $oldstddev = sprintf("%s.%s", $5, $6) + 0.0;
        if (exists $ResultLibrary{$acc})
        {
          $ResultLibrary{$acc}{lid}  = $acc;         # cross reference
          $ResultLibrary{$acc}{name}  = $name;         # cross reference
          $InputLibrary{$acc}{name}  = $name;
          $InputLibrary{$acc}{mean} = $oldmean;
          $InputLibrary{$acc}{stddev}= $oldstddev;
          $ResultLibrary{$acc}{stats}{input_mean} = $oldmean;
          $ResultLibrary{$acc}{stats}{input_stddev}= $oldstddev;
        }
        else
        {
          logwarning("Did not encounter library $name (lid=$acc) in .asm file.");
        }
      }

      # Emit rewritten .insert record if specified
      ++$nrecords;
      progress("Phase 1: $nrecords records scanned.") if ($nrecords > 1 && $nrecords % $PROGRESS_FACTOR == 1);
  }

  $if->close();

  # Apply filters
  # Test 1: Remove any libraries for which there are no experimental data
  #
  if ($options{removeEmpty})
  {
    my @empties = ();
    foreach my $acc (keys %ResultLibrary)
    {
      my $N = (exists $ResultLibrary{$acc}{stats}{items})? $ResultLibrary{$acc}{stats}{items} : 0;
      my $name = (exists $ResultLibrary{$acc}{name})? $ResultLibrary{$acc}{name} : "name?";
      if ($N == 0)
      {
        push(@empties, $acc) if ($N == 0);
        logwarning("Library $acc (\'$name\') has no results, omitting from output.");
      }
    }
    foreach my $acc (@empties)
    {
      delete $ResultLibrary{$acc};
    }
  }

  # Setup output streams
  if (defined $options{cloneoutfile})
  {
    progress("Writing .clone file \'$options{cloneoutfile}...");
    my $of = new IO::File("> $options{cloneoutfile}")
      or $tf->bail("Cannot open output .clone file \'$options{cloneoutfile}\' ($!)");
    foreach my $acc (keys %ResultLibrary)
    {
      my $name = (exists $InputLibrary{$acc}{name})? $InputLibrary{$acc}{name} : "name?";
      my $mean = $ResultLibrary{$acc}{stats}{mean};
      my $stddev = $ResultLibrary{$acc}{stats}{stddev};
      if (exists $InputLibrary{$acc}{name})
      {
        $of->print("$name\t$mean $stddev\n");
      }
      else
      {
        logwarning("Unknown library encountered in .asm file: $name (lid=$acc, mea=$mean, std=$stddev). Omitting from .clone file...");
      }
    }
    $of->close() or $tf->bail("Failed to close output .clone file ($!)");
    progress("Done.");
  }

  if (defined $options{insertoutfile})
  {
    progress("Writing new .insert file \'$options{insertoutfile}\'...");
    my $if = new IO::File("< $options{insertfile}")
      or $tf->bail("Cannot open input .insert file \'$options{insertfile}\' ($!)");
    my $sf = new IO::File("> $options{insertoutfile}")
      or $tf->bail("Cannot open output .insert file \'$options{insertoutfile}\' ($!)");

    my $nrecords = 0;
    while (my $line = $if->getline())
    {
      chomp $line;
      # I wish we had a proper module to interact with this file type
      if ($line =~ m/\s+<GenomicLibrary Id="(\w+)" acc="(\d+)" mea="(\d+).(\d+)" std_dev="(\d+).(\d+)">/o)
      {
        my $name = $1;
        my $acc = $2;
        if (exists $ResultLibrary{$acc})
        {
          my $line2 = "";
          $line2 .= "  <GenomicLibrary ";
          $line2 .= "Id=\"$ResultLibrary{$acc}{name}\" ";
          $line2 .= "acc=\"$acc\" ";
          $line2 .= "mea=\"$ResultLibrary{$acc}{stats}{mean}\" ";
          $line2 .= "std_dev=\"$ResultLibrary{$acc}{stats}{stddev}\">";
          $sf->print("$line2\n") or $tf->bail("Cannot write output .insert file \'$options{insertoutfile}\' ($!)");
        }
        else
        {
          # pass through
          $sf->print("$line\n") or $tf->bail("Cannot write output .insert file \'$options{insertoutfile}\' ($!)");
        }
      }
      else
      {
        $sf->print("$line\n") or $tf->bail("Cannot write output .insert file \'$options{insertoutfile}\' ($!)");
      }

      # Emit rewritten .insert record if specified
      ++$nrecords;
      progress("Phase 1: $nrecords records scanned.") if ($nrecords > 1 && $nrecords % $PROGRESS_FACTOR == 1);
    }
    $if->close();
    $sf->print("\n");   # make sure last record is set off
    $sf->close() or $tf->bail("Failed to close output .insert file ($!)");
    progress("Done.");
  }

  if ($options{metricsfile})
  {
    progress("Writing .metrics file \'$options{metricsfile}\'...");
    my $mf = new IO::File("> $options{metricsfile}")
      or $tf->bail("Cannot open output .metrics file \'$options{metricsfile}\' ($!)");

    foreach my $acc (keys %ResultLibrary)
    {
      my $line = "";
      my $name = $ResultLibrary{$acc}{name};
      $line .= "[library_$acc]\n";                     # title
      $line .= "lid=$acc\n";
      $line .= "name=$ResultLibrary{$acc}{name}\n";
      if (exists $ResultLibrary{$acc}{stats}{input_mean})
      {
        $line .= "input_mean=$ResultLibrary{$acc}{stats}{input_mean}\n";
      }
      else
      {
        logwarning("Did not find existing mean for library $name (lid=$acc).");
      }
      $line .= "mean=$ResultLibrary{$acc}{stats}{mean}\n";
      if (exists $ResultLibrary{$acc}{stats}{input_stddev})
      {
        $line .= "input_sd=$ResultLibrary{$acc}{stats}{stddev}\n";
      }
      else
      {
        logwarning("Did not find existing stddev for library $name (lid=$acc).");
      }
      $line .= "sd=$ResultLibrary{$acc}{stats}{stddev}\n";
      $line .= "min=$ResultLibrary{$acc}{stats}{min}\n";
      $line .= "max=$ResultLibrary{$acc}{stats}{max}\n";
      $line .= "N=$ResultLibrary{$acc}{stats}{items}\n";
      $line .= "binsize=$ResultLibrary{$acc}{hist}{bin}\n";
      $line .= "buc=$ResultLibrary{$acc}{hist}{buc}\n";
      $line .= "hist=" . join(",",@{$ResultLibrary{$acc}{hist}{values}}) . "\n";
      $mf->print("$line\n") or $tf->bail("Failed to write output .metrics file \'$options{metricsfile}\' ($!)");
    }
    $mf->close() or $tf->bail("Failed to close output .metrics file ($!)");
  }

  if ($options{libraryoutfile})
  {
    progress("Writing .library file \'$options{libraryoutfile}\'...");
    my $xf = new IO::File("> $options{libraryoutfile}")
      or $tf->bail("Cannot open output .library file \'$options{libraryoutfile}\' ($!)");

    my $header        = qq~<?xml version='1.0'?>~ . "\n";
    my $group         = "<librarystats>\n";
    my $groupend      = "</librarystats>\n";
    my $library       = "\t<library>\n";
    my $libraryend    = "\t</library>\n";
    my $stats         = "\t\t<stats>\n";
    my $statsend      = "\t\t</stats>\n";
    my $histogram     = "\t\t<histogram>\n";
    my $histogramend  = "\t\t</histogram>\n";
    my $bucketdataF   = qq~buckets="%d"~;
    my $bucketsizeF   = qq~bucketsize="%0.2f"~;

    my $sectionF      = "\t\t<%s>%s</%s>\n";
    my $histogramF    = "\t\t<histogram>\n\t\t\t<listdata buckets=\"%s\" bucketsize=\"%s\">\n\t\t\t\t%s\n\t\t\t</listdata>\n\t\t</histogram>\n";
    my $statslineF    = "\t\t\t<%s>%s</%s>\n";

    $xf->print ($header) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
    $xf->print ($group) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
    foreach my $acc (keys %ResultLibrary)
    {
      $xf->print ($library) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");

      my $rh_value2section = rankOrder(\%{$ResultLibrary{$acc}}, \%SECTION_ITEMS);
      foreach my $rank (sort { $a<=>$b } keys %{$rh_value2section})
      {
        my $section = ${$rh_value2section}{$rank};
        if ($section eq "stats")
        {
          $xf->print ($stats) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
          {
            my $rh_value2stat = rankOrder(\%{$ResultLibrary{$acc}{stats}}, \%STAT_ITEMS);
            foreach my $rank (sort { $a<=>$b } keys %{$rh_value2stat})
            {
              my $tag = ${$rh_value2stat}{$rank};
              my $line = sprintf($statslineF, $tag, $ResultLibrary{$acc}{stats}{$tag}, $tag);
              $xf->print ($line) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
            }
          }
          $xf->print ($statsend) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
        }
        elsif ($section eq "hist")
        {
          my $line = sprintf($histogramF, $ResultLibrary{$acc}{hist}{buc}, $ResultLibrary{$acc}{hist}{bin},
                                          join(",", @{$ResultLibrary{$acc}{hist}{values}}) );
          $xf->print ($line) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
        }
        else
        {
          my $line = sprintf($sectionF, $section, $ResultLibrary{$acc}{$section}, $section);
          $xf->print ($line) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
        }
      }
      $xf->print ($libraryend) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
    }
    $xf->print ($groupend) or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");
    $xf->print ("\n") or $tf->bail("Failed to write output .library file \'$options{libraryoutfile}\' ($!)");

    $xf->close() or $tf->bail("Failed to close output .library file ($!)");
  }

  progress("Done.");
  exit 0;
}

