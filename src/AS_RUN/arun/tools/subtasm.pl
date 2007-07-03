#!/usr/local/bin/perl
#
# subtasm - Filter, subset, or merge .tasm records depending on input
#           parameters.
#
# Written by: Martin Shumway
#

use IO::File;
use File::Basename;
use TIGR::Foundation;
use TIGR::AsmLib;
use TIGR::EUIDService;
use strict;
use warnings;

my $tf = new TIGR::Foundation;
die "Could not establish TIGR::Foundation object" if (! defined $tf);
setFoundation($tf);
my $te = new TIGR::EUIDService;
die "Could not establish TIGR::EUIDService object" if (! defined $te);

my $PRG = $tf->getProgramInfo('name');
my $VERSION="1.12";
my @DEPENDS=("TIGR::Foundation");

my $HELPTEXT = qq~
Extract, filter, and merge .tasm files from an CA or TA assembly run. 

subtasm  <tasm1 tasm2 ...>  [options]
    
  tasms  One or more .tasm files as constructed by TA/CA pipeline and analysis
         tools.  File list can be resolved by the shell, for example my*.tasm
    
  options:
    -A <file>      File with cids (ca_contig_id EUIDs) to extract for output
    -X <file>      File with cids (ca_contig_id EUIDs) to omit from output
    -x <file>      Strip these reads from their contigs and recompute seq#
    -r <file>      Find and replace according to perl regexp rule.  Regexp spec
                   file is tab delimited with the first field as the tag name
                   and second field as the search/replace rule.
    -f <file>      Substitute existing seq_name (eg sequence read EUID) 
                   according to the input lookup file.  Format is <rid><alias> 
                   as is encoded in the .seq.features file produced by pullfrag.
    -o <files>     Write output to file instead of the console.
    -[no]circular  Set circular field for contigs
    -comment <s>   Set comment field for contigs
    -com_name <s>  Set com_name field for contigs
 
Given one or more input .tasm files (or .asm files produced by run_TA) 
subtasm builds a single new .tasm file with suitable treatments.

If no filter option is specified then all records are emitted on output.  
Specifying multiple files on the command line has the effect of merging the
tasm records prior to filtering.  The output contig set can be selected (-A) or
decimated (-X).  On output field values can be modified according to a lookup
table (-f) or a regular expression substitution (-r). Finally, individual reads
can be removed from their placements (-x) recomputing seq#.  If the result
is 0, these will be dropped from the output file.

subtasm is useful for performing in-situ operations on .tasm files prior to
uploading to the database.  Additional output conversions can be specified
using the -o option.  Otherwise, a new .tasm file is generated on stdout.

SEE ALSO
  aloader carun run_TA slice2tasm 
~;

# =============================== Constants ================================

# Operational parameters
my $debug = 0;
my $LOW = 1;
my $MEDIUM = 2;
my $HIGH = 3;
my $PROGRESS_FACTOR = 100000;
my $warningsfile = "tmp.warnings";
my $USERNAME = getpwuid($<);

my @TASM_CONTIG_TAGS = 
(
  "sequence",
  "lsequence",
  "asmbl_id",
  "seq_id",
  "com_name",
  "type",
  "method",
  "ed_status",
  "redundancy",
  "perc_N",
  "seq#",
  "full_cds",
  "cds_start",
  "cds_end",
  "ed_pn",
  "ed_date",
  "comment",
  "frameshift",
  "ca_contig_id",
  "mod_date",
  "mod_pn",
  "is_circular",
);

my @TASM_READ_TAGS = 
(
  "seq_name",
  "asm_lend",
  "asm_rend",
  "seq_lend",
  "seq_rend",
  "best",
  "comment",
  "db",
  "offset",
  "lsequence",
);

# Supports use of .asm files from TIGR Assembler (run_TA) as well as standard
# timmy files (.tasm)
my @OUTPUT_SUFFIXES = ("asm", "tasm");  
my %SUPPORTED_OUTPUTS = ();
map { $SUPPORTED_OUTPUTS{$_} = 1; } @OUTPUT_SUFFIXES;

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

sub emitContigRecord($$$$$$$$$$$)
{
    $tf->logLocal("emitContigRecord");
  my ($of, $rh_CurrentContig, $rh_options, $rh_FilterSet, $current_contig, $first, $rh_ContigsStripped, $rh_Contigs, $comment, $circular, $com_name) = @_;
  my $emit = 1;

  if (${$rh_options}{filter})
  {
    if (${$rh_options}{FilterByIncludeList})
    {
      $emit = (exists ${$rh_FilterSet}{$current_contig})? 1 : 0;
    }
    if (${$rh_options}{FilterByExcludeList})
    {
      $emit = (exists ${$rh_FilterSet}{$current_contig})? 0 : 1;
    }
  }
  else
  {
    $emit = 1;   # Pass through
  }

  # Suppress if the contig has no reads in it.
  if ( exists ${$rh_ContigsStripped}{$current_contig} )
  {
    $emit = 0;
  }

  return if (! $emit);

  unless ($first)
  {
    $of->print("|\t\n") or $tf->bail("Failed to write output tasm stream ($!)");
  }

  # fix the seq count to what was discovered
  ${$rh_CurrentContig}{'seq#'} = ${$rh_Contigs}{$current_contig};

  my ($sec,$min,$hour,$mday,$mon,$year)=localtime(time);

  my $amPm = ($hour >= 12 and $hour < 24) ? 'PM' : 'AM';
  $hour %= 12;
  $hour = 12 if ( $hour == 0 );

  my $dateString = sprintf("%02d/%02d/%02d %02d:%02d:%02d %s",$mon+1,$mday,$year%100,$hour,$min,$sec,$amPm);
   
  foreach my $t (@TASM_CONTIG_TAGS)
  {
    my $printStr;
    if ($t eq "comment" && defined $comment ) {
      ${$rh_CurrentContig}{$t} = $comment;
      $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
    } elsif ( $t eq "com_name" && defined $com_name ) {
      ${$rh_CurrentContig}{$t} = $com_name;
      $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
    } elsif ( $t eq "is_circular" && defined $circular ) {
      ${$rh_CurrentContig}{$t} = $circular;
      $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
    } elsif (exists ${$rh_CurrentContig}{$t} ) {
      if ($t eq "ca_contig_id") 
      {
        if (! defined ${$rh_CurrentContig}{$t}  ||  
            (defined ${$rh_CurrentContig}{$t} && ${$rh_CurrentContig}{$t} eq "")
           )
        {
          # ca_contig_id (cid) exists, but is empty    
          my $euid = $te->getEUID();
          ${$rh_CurrentContig}{$t} = $euid;
        }
      }
      elsif ($t eq "ed_pn"  &&  
           ( 
             (! defined ${$rh_CurrentContig}{$t})  ||  
             (defined ${$rh_CurrentContig}{$t} && ${$rh_CurrentContig}{$t} eq "GRA") 
           )
         )
      {
        ${$rh_CurrentContig}{$t} = $USERNAME; 
      }
      elsif ($t eq "is_circular" && ${$rh_CurrentContig}{$t} =~ /^\s*$/ ) 
      {
          ${$rh_CurrentContig}{$t} = 0;
      }
      elsif ($t eq "comment" && ${$rh_CurrentContig}{$t} =~ /^CA_FREE/ ) 
      {
          ${$rh_CurrentContig}{$t} =~ s/CA_FREE/Non redundified surrogate contig\. CA_FREE/;
      }
      elsif ($t eq "method" && ${$rh_CurrentContig}{$t} =~ /Celera Assembler/ ) 
      {
          ${$rh_CurrentContig}{$t} = "CA ";
      }
      elsif ($t eq "mod_date" ) {
          ${$rh_CurrentContig}{$t} = $dateString;
          $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
      }
      $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
    }
    elsif ($t eq "ca_contig_id") 
    {
      # ca_contig_id (cid) does not exist so we have to assign one
      my $euid = $te->getEUID();
      ${$rh_CurrentContig}{$t} = $euid;
      $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
    }
    elsif ($t eq "is_circular" ) {
      ${$rh_CurrentContig}{$t} = 0;
      $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
    }
    elsif ($t eq "mod_pn" ) {
      ${$rh_CurrentContig}{$t} = $USERNAME;
      $printStr="$t\t${$rh_CurrentContig}{$t}\n";      
    }
    else
    {
      $printStr="$t\t\n";      
    }
    $of->print($printStr) or $tf->bail("Failed to write output tasm stream ($!)");
  }
}

sub emitReadRecord($$$$$$$)
{
  my ($of, $rh_CurrentRead, $rh_options, $rh_FilterSet, $current_contig, $current_read, $rh_StripSet) = @_;
  my $emit = 1;

  # Contig filtering by inclusion or exclusion
  if (${$rh_options}{filter})
  {
    if (${$rh_options}{FilterByIncludeList})
    {
      $emit = (exists ${$rh_FilterSet}{$current_contig})? 1 : 0;
    }
    if (${$rh_options}{FilterByExcludeList})
    {
      $emit = (exists ${$rh_FilterSet}{$current_contig})? 0 : 1;
    }
  }
  else
  {
    $emit = 1;   # pass through
  }

  # Read filtering by exclusion
  if (${$rh_options}{strip})
  {
    $emit = (exists ${$rh_StripSet}{$current_read})? 0 : $emit;
  }

  return if (! $emit);

  $of->print("\n") or $tf->bail("Failed to write output tasm stream ($!)");
  foreach my $t (@TASM_READ_TAGS)
  {
    if (exists ${$rh_CurrentRead}{$t})
    {
      $of->print("$t\t${$rh_CurrentRead}{$t}\n") or $tf->bail("Failed to write output tasm stream ($!)");
    }
    else
    {
      $of->print("$t\t" . "") or $tf->bail("Failed to write output tasm stream ($!)");
    }
  }
}

# ============================================== MAIN =============================================
#
MAIN:
{    
  my %options = ();
  $options{FilterByIncludeList}   = undef;
  $options{FilterByExcludeList}   = undef;
  $options{FilterByStripList}     = undef;
  $options{FilterByFindReplace}   = undef;
  $options{FilterBySubstitution}  = undef;
  $options{outprefix}             = undef;
  $options{comment}               = undef;
  $options{circular}              = undef;
  $options{com_name}              = undef;
  $options{filter}     = 0;
  $options{filterfile} = undef;
  $options{substitute}     = 0;
  $options{substitutionfile} = undef;
  $options{replace}     = 0;
  $options{replacefile} = undef;

  # These tables need to be populated regardless of input filtering modes.
  # 
  my %Catalog         = (); # Contig to read lookup
  my %FilterSet       = (); # List of contigs to extract
  my %StripSet        = (); # List of reads to strip    
  my %RegexpRules     = (); # List of tags and perl regexp rules 
  my %SubstitutionSet = (); # List of objects to substitute
  my %Contigs         = (); # List of contigs and their seq# settings
  my %ContigsStripped = (); # List of contigs with no reads
  
  # Configure TIGR Foundation
  $tf->setHelpInfo($HELPTEXT);
  $tf->setUsageInfo($HELPTEXT);
  $tf->setVersionInfo($VERSION);
  $tf->addDependInfo(@DEPENDS);
  
  # validate input parameters
  my $output_options = undef;
  my $result = $tf->TIGR_GetOptions
               (
                'A:s'       => \$options{FilterByIncludeList}, 
                'X:s'       => \$options{FilterByExcludeList},
                'x:s'       => \$options{FilterByStripList},
                'r:s'       => \$options{FilterByFindReplace},                         
                'f:s'       => \$options{FilterBySubstitution},                         
                'o:s'       => \$options{outprefix},
                'comment=s' => \$options{comment},
                'circular!' => \$options{circular},
                'com_name=s'=> \$options{com_name},
               );
               print "Result=$result\n";
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

  my $warningsfile = (defined $options{outprefix})?  "$options{outprefix}.warnings" : "$PRG.warnings";
  unlink $warningsfile;
  setWarnFile($warningsfile); 

  unless ($te->ping()) 
  {
    $tf->bail("EUID service failed ping test and is therefore unavailable, exiting...");
  }

  my $ncontigs = 0;
  my @infiles = ();
  my $nfiles = $#ARGV + 1;
  for (my $i=0; $i <= $#ARGV; $i++)
  {
    my $infilename =  $ARGV[$i];
    push @infiles, $infilename; 
    $tf->bail("Cannot access input file \'$infilename\' ($!)") if (! -r $infilename);
    my $n_asmbl_id= `grep -c asmbl_id $infilename`; 
    $ncontigs += $n_asmbl_id;
  }
  progress("Encountered $ncontigs contigs in $nfiles files.");

  if (defined $options{FilterByIncludeList})
  {
    $options{filterfile} = $options{FilterByIncludeList};
    $options{filter} = 1;
  }
  if (defined $options{FilterByExcludeList})
  {
    $options{filterfile} = $options{FilterByExcludeList};
    $options{filter} = 1;
  }
  if (defined $options{FilterByStripList})
  {
    $options{stripfile} = $options{FilterByStripList};
    $options{strip} = 1;
  }
  if (defined $options{FilterByFindReplace})
  {
    $options{replacefile} = $options{FilterByFindReplace};
    $options{replace} = 1;
  }
  if (defined $options{FilterBySubstitution})
  {
    $options{substitutionfile} = $options{FilterBySubstitution};
    $options{substitute} = 1;
  }

  # Get optional list of filter of sequences to output.  If this
  # option is not supplied then the whole input set is used.
  #
  if ($options{filter})
  {
    #get read names from input file
    progress("Obtaining filter rules from \'$options{filterfile}\'");
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
  else
  {
    progress("No filter specified, using all contigs.");
  }

  if ($options{strip})
  {
    #get read names from input file
    progress("Obtaining strip list from \'$options{stripfile}\'");
    my $ff = new IO::File("< $options{stripfile}")
      or $tf->bail("Failed to open input strip file \'$options{stripfile}\' ($!)");
    while (my $line = $ff->getline())
    {
      chop $line;
      next if ($line =~/^#/);
      my @f = split /\t/,$line;
      $StripSet{$f[0]} = $f[1]; 
    }
    $ff->close();
  }
  else
  {
    progress("No reads to strip, using all reads.");
  }

  if ($options{replace})
  {
    #get read names from input file
    progress("Obtaining regexp rules from \'$options{replacefile}\'");
    my $ff = new IO::File("< $options{replacefile}")
      or $tf->bail("Failed to open input regexp rules file \'$options{replacefile}\' ($!)");
    while (my $line = $ff->getline())
    {
      chop $line;
      next if ($line =~/^#/);
      my @f = split /\t/,$line;
      my $tag = $f[0];
      my $value = (defined $f[1])? $f[1] : "";
      $RegexpRules{$tag} = $value;
    }
    $ff->close();
  }
  else
  {
    progress("No regexp find/replace specified, modifying no records.");
  }

  if ($options{substitute})
  {
    #get read names from input file
    progress("Obtaining substitution rules from \'$options{substitutionfile}\'");
    my $ff = new IO::File("< $options{substitutionfile}")
      or $tf->bail("Failed to open input filter file \'$options{substitutionfile}\' ($!)");
    while (my $line = $ff->getline())
    {
      chop $line;
      next if ($line =~/^#/);
      my @f = split /\t/,$line;
      $SubstitutionSet{$f[0]} = $f[1]; 
    }
    $ff->close();
  }
  else
  {
    progress("No substition specified, modifying no records.");
  }

  progress("Screen input files to build contig-read catalog...");
  my $current_contig_id = undef;
  foreach my $infile (@infiles)
  {
    my $if = new IO::File("< $infile") or $tf->bail("Cannot open input .tasm file \'$infile\' ($!)");
    while (my $line = $if->getline())
    {
      chop $line;
      my @f = split /\t/,$line;
      next if (! defined $f[0]);
      if ($f[0] eq "asmbl_id"  ||  $f[0] eq "ca_contig_id")
      {
        if (defined $f[1]  &&  $f[1] ne "")
        {
          $current_contig_id = $f[1];
          $Contigs{$current_contig_id} = 0;  # initialize read count
        }
      }
      elsif ($f[0] eq "seq_name")
      {
        $Catalog{$f[1]} = $current_contig_id;
        $Contigs{$current_contig_id}++;
      } 
    } 
    $if->close(); 
  }

  foreach my $xread (keys %StripSet)
  {
    if (exists $Catalog{$xread})
    {
      my $cid = $Catalog{$xread};
      $Contigs{$cid}--;
      delete $Catalog{$xread};
    }
    else
    {
      logwarning("Strip read \'$xread\' not found in contig set, ignoring.");
    }
  } 

  foreach my $contig (keys %Contigs)
  {
    if ($Contigs{$contig} <= 0)
    {
      $ContigsStripped{$contig} = 1;
      delete $Contigs{$contig};
      logwarning("Contig \'$contig\' has no reads, removing.");
    }
  }

  # Setup output streams
  $options{tasmout} = (defined $options{outprefix})? "$options{outprefix}.tasm" : undef;
  my $of = undef;
  if (defined $options{tasmout})
  {
    my $outfilename = $options{tasmout};
    $of = new IO::File("> $outfilename") or $tf->bail("Cannot open output file \'$outfilename\' ($!)");
    progress("Output being sent to \'$outfilename\'");
  }
  else
  {
    $of = new IO::Handle;
    $of->fdopen(fileno(STDOUT), "w");
    progress("Output being sent to stdout"); 
  }

  # Stream out the input file such that each FRG record is read and emitted
  # with possible modifications depending on execution mode.
  #
  progress("Process and write output...");
  my %CurrentContig = ();    
  my $current_contig = undef;
  my %CurrentRead = ();      
  my $current_read = undef;
  my $inContig = 1;
  my $nseqs = 0;
  my $first = 1;
  my $if = undef;
  foreach my $infile (@infiles)
  {
    if (defined $if  &&  $if->eof())
    {
      # Reset to be ready to take on next contig record
      %CurrentRead = ();    # flush
      $current_read = undef;
      %CurrentContig = ();  # flush
      $current_contig = undef;
      $inContig = 1;        # next line is in a contig record
      $nseqs = 0; 
    }
    $if = new IO::File("< $infile") or $tf->bail("Cannot open input .tasm file \'$infile\' ($!)");
    progress("Processing input file \'$infile\'...");
    my $nline = 0;
    while (my $line = $if->getline())
    {
      chop $line;
      $nline++;
      my @f = split /\t/,$line;
      my $tag = $f[0];
      my $value = (defined $f[1])? $f[1] : "";

      if (defined $tag)
      {     
        # Start of new contig, flush existing read record
        if ($tag eq "|") 
        {
          # What if a leading or trailing contig delimiter found
          if (! defined $current_contig)
          {
            next;  # eat leading or trailing delimiters in file
          }
          emitReadRecord($of , \%CurrentRead, \%options, \%FilterSet, $current_contig, $current_read, \%StripSet);
          %CurrentRead = ();
          $current_read = undef;
          $inContig = 1;        # next line is in a contig record
          $nseqs++; 
          next;
        }

        # All other tag-value pairs
        else
        {
          # else bind tag value pairs with treatments
          # Bind current record(s)
          my $current_contig1 = ($tag eq "asmbl_id"     && defined $value && $value ne "")? $value : undef;
          my $current_contig2 = ($tag eq "ca_contig_id" && defined $value && $value ne "")? $value : undef;
          $current_contig = (defined $current_contig1)? $current_contig1 : 
                            (defined $current_contig2)? $current_contig2 : 
                            $current_contig;
          $current_read = ($tag eq "seq_name" && defined $value && $value ne "")? $value : $current_read;
    
          # Apply find/replace rule to each record
          if (exists $RegexpRules{$tag}  &&  $options{FilterByFindReplace})
          {
            my $expr = "\$value =~ " . $RegexpRules{$tag} . ";";
            my $result = eval($expr); 
            if (! $result)
            {
              $tf->logLocal("Failed to apply expression \'$expr\' to field \'$tag\': ($@): current value = $value", 3, 1);
            }
          }

          # Apply substitution rule to every seq_name record
          if ($tag eq "seq_name"  &&  $options{FilterBySubstitution})
          {
            $value = $SubstitutionSet{$value}; 
          }

          if ($inContig)
          {
            $CurrentContig{$tag} = (defined $value)? $value : "";
          }
          else
          {
            $CurrentRead{$tag} = (defined $value)? $value : "";
          }

          # Emit residual read record
          if ($if->eof())
          {
            # Spit out last record in queue following end of file
            emitReadRecord($of , \%CurrentRead, \%options, \%FilterSet, $current_contig, $current_read, \%StripSet);
            %CurrentRead = ();
            $inContig = 0;        # next line is in a read record
            $nseqs++; 
          }
        }
      }
      
      # Tag is not defined.  Start readset or inter-read delimiter 
      else      
      {
        # Another read record, flush current read record
        if (defined $current_read)
        {
          if (exists $Catalog{$current_read} )
          {
            emitReadRecord($of , \%CurrentRead, \%options, \%FilterSet, $current_contig, $current_read, \%StripSet);
          }
          %CurrentRead = ();
          $inContig = 0;        # next line is in a read record
          $nseqs++; 
          next;
        }
   
        # Read set about to begin, flush current contig record
        else
        {
          emitContigRecord($of, \%CurrentContig, \%options, \%FilterSet, $current_contig, $first, \%ContigsStripped, \%Contigs, $options{comment}, $options{circular}, $options{com_name});
          %CurrentContig = ();
          $inContig = 0;        # next line is in a read record
          $ncontigs++;
          $first = 0;
          next;
        }
      }
    }

    $if->close();
  }

  # Close output stream
  if (defined $of)
  {
    $of->close() or $tf->bail("Failed to close output frg stream ($!)");
  }

  progress("Done.");
  exit 0;
}

