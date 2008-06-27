#!/usr/local/bin/perl -w
# $Id: tarun.pl,v 1.3 2008-06-27 06:29:19 brianwalenz Exp $
#
# TIGR Assembler front-end script calls/runs assembly pipline
# recipes on the grid or on the local host.
#
# Written by Prabhu May 2007.
#

# Program configuration
# ================================== Constants  ================================
# Program configuration
my $MY_APP = basename($0);
my $MY_LOG = "$MY_APP.log";
my $MY_VERSION = " Version 1.00";

my $MY_DATE_FORMAT = "%4d_%02d%02d";
my $MY_APPLICATION = "TIGR_Assembler";
#my $MLIMIT_FACTOR = 72;
my $DEFAULT_MIN_PCT = 95;
##my $RUNCMD_EXEC = "/home/sgeworker/bin/runCmd";
my $RUNCMD_EXEC = "runCmd";
my $NOTIFY_SCRIPT = "/usr/local/common/tarun_notify";
my $CMD_LOG = "tarun.log";

my @MY_DEPENDS = ( "TIGR::Foundation", "TIGR::ConfigFile" );
my $HELPTEXT =
  qq~Request a whole-genome shotgun assembly using TIGR Assembler.

 usage: tarun [options] <prefix>
	    tarun -cancel <id>

 inputs:
  <prefix>       Prefix for sequence and qual files

 general options:
  -alias <a>        Identify the request as <a> on Assembly Server Console
  -cancel <id> 	    Cancel assembly request <id>
  -maxCopy <dir>    Copy max output to <dir> (default: dir=maxCopy)
  -medCopy <dir>    Copy med output to <dir> (default: dir=medCopy)
  -minCopy <dir>    Copy min output to <dir> (default: dir=minCopy)
  -noCopy           Do not make local copy. (default: medCopy)
  -local            Combination of 'localHost' and 'localDisk'
  -localHost        Assemble on this host (default: on grid)
  -localDisk        Use current directory to assemble. (default: aserver)
  -D <name>         Specify project name (displayed on AserverConsole).
                    Use -D test for debugging. (default: frg prefix)
  -[no]notify    	Send email upon completion/error and create BITS case
                    for errors. (default: notify)
  -wait             Wait for job completion before exiting

 recipe options:
  -R <standard_recipe to execute TIGR Assembler> OR TA2 (default) OR TA2.5

 tarun [-local] [options] file
    -local  Force the TIGR Assembler to run on the local system with output
            going to the file <file>.tasm .  The default behavior is to submit
            the assembly job for execution on a system allocated from the
            sge pool, and receive email notification that the job is done.
            Use the -local option for small jobs only.

  	-recipe <recipe>       Specify recipe to execute TIGR Assembler

    -dir    Specify another directory pathname for the output files
 TA   options
        -q <qual file>  quality file or pathname
        -C <contig file>  Contig jump start or pathname
        -[no]constraints  Apply constraints to assembly process.
        -e <n>  Set the max length that a DNA fragment may overhang
                to be considered for assembly.
        -l <n(int)>  Set the min length by which two DNA fragments must
                overlap to be considered for assembly.
        -p <f(float)>  Set the min percent identity that two DNA fragments must
                achieve in their overlap region to be considered for assembly.
        -s      Include singletons (single sequence contigs) in output.
        -bigjob   For jobs exceeding grid's fast queue limit.
        -o		option for output
    file     Name of the multifasta formatted .seq input file of DNA sequences.


 Tracking assembly requests:
 http://assemblyconsole.tigr.org/
~;
my $HELPTEXT_BETA = qq~
Request a whole-genome shotgun assembly using TIGR Assembler.

 usage: tarunb [options] <prefix>
	    tarunb -cancel <id>

 inputs:
  <prefix>       Prefix for sequence and qual files

 general options:
  -alias <a>        Identify the request as <a> on Assembly Server Console
  -cancel <id> 	    Cancel assembly request <id>
  -maxCopy <dir>    Copy max output to <dir> (default: dir=maxCopy)
  -medCopy <dir>    Copy med output to <dir> (default: dir=medCopy)
  -minCopy <dir>    Copy min output to <dir> (default: dir=minCopy)
  -noCopy           Do not make local copy. (default: medCopy)
  -local            Combination of 'localHost' and 'localDisk'
  -localHost        Assemble on this host (default: on grid)
  -localDisk        Use current directory to assemble. (default: aserver)
  -D <name>         Specify project name (displayed on AserverConsole).
                    Use -D test for debugging. (default: frg prefix)
  -[no]notify    	Send email upon completion/error and create BITS case
                    for errors. (default: notify)
  -wait             Wait for job completion before exiting

 recipe options:
  -R <standard_recipe to execute TIGR Assembler> OR TA2 (default) OR TA2.5

 tarun [-local] [options] file
    -local  Force the TIGR Assembler to run on the local system with output
            going to the file <file>.tasm .  The default behavior is to submit
            the assembly job for execution on a system allocated from the
            sge pool, and receive email notification that the job is done.
            Use the -local option for small jobs only.

  	-recipe <recipe>       Specify recipe to execute TIGR Assembler

    -dir    Specify another directory pathname for the output files
 TA2.5   options
        -q <qual file>  quality file or pathname
        -C <contig file>  Contig jump start or pathname
        -[no]constraints  Apply constraints to assembly process.
        -e <n>  Set the max length that a DNA fragment may overhang
                to be considered for assembly.
        -l <n(int)>  Set the min length by which two DNA fragments must
                overlap to be considered for assembly.
        -p <f(float)>  Set the min percent identity that two DNA fragments must
                achieve in their overlap region to be considered for assembly.
        -s      Include singletons (single sequence contigs) in output.
        -U <f>  Set the start percent for stepwise assembly (-U -B -i combo).
        -B <f>  Set the end percent for stepwise assembly (-U -B -i combo).
        -i <f>  Set the increment percent for stepwise assembly (-U -B -i combo).

        -o		option for output
    file     Name of the multifasta formatted .seq input file of DNA sequences.


 Tracking assembly requests:
 http://assemblyconsole.tigr.org/
 ~;


# ================================== Pragmas ==================================
use strict;
use POSIX qw(uname);
use POSIX qw(:signal_h :errno_h :sys_wait_h);
use sigtrap 'handler', \&SIGHANDLER, 'normal-signals', 'error-signals',
  'old-interface-signals';
use TIGR::Foundation;
use TIGR::ConfigFile;
use Getopt::Long;
use IO::File;
use File::Basename;
# ============================== Program Variants ==============================

# Which TIGR_Assembler image to run.

my %APPLICATIONS =
(
  "tarun" =>      "TIGR_Assembler",
  "tarun.pl" =>   "TIGR_Assembler",
  "tarun_beta" => "TIGR_Assembler_beta_ifo",
);

my $TIGRASM         = "/usr/local/packages/$MY_APPLICATION/bin/$APPLICATIONS{$MY_APP}";
my $LOCAL_TIGRASM   = "$TIGRASM";

# What option string to give TIGR_Assembler
my %OPTIONS =
(
  "tarun" =>      "-g 8 -l 40 -e 15 -p 97.5",                # These are overridden by command line
  "tarun_beta" =>  "-g 8 -l 40 -e 15 -p 97.5 -constraints",   # These are overridden by command line
);
my $MY_OPTIONS = $OPTIONS{$MY_APP};
# Does TIGR_Assembler support safe resource allocation?
my %SAFEOPS =
(
  "tarun" =>      0,
  "tarun_beta" =>  1,     # add additional options to ensure safe resource allocation
);
# Does TIGR_Assembler give its version (supports -V option)?
my %HAS_VERSION =
(
  "tarun" =>      0,
  "tarun_beta" =>  1,     # has -V option
);

my $MY_HASVERSION = $HAS_VERSION{$MY_APP};

my $APP_OPTIONS = "";
# Process the application-specific options
# parm 0 - The option
# parm 1 - Additional option, if supplied
sub appOptions
{
  $APP_OPTIONS .= " -$_[0]";
  $APP_OPTIONS .= " $_[1]" if (defined $_[1]);
}

# Process the application-specific flag
# parm 0 - The option
sub appFlags
{
  $APP_OPTIONS .= " -$_[0]";
}

my $MY_SAFEOPS = $SAFEOPS{$MY_APP};
# ================================== Constants  ===============================
# default config file for tarun
my $CARUN_CONFIG_FILE = "/usr/local/common/ARUN/CARUN.conf";

# The tarun Config file object
my $carun_cf = undef;

# when Cntl-C is used to cancel tarun
my $cancel_track = 0;

# Module for housekeeping functions
my $tf_object = new TIGR::Foundation;

# Stack of files to clean at the end of the run
my @filesToClean = ();


#Name:   setRecipeFile
#Input:  recipe_file per user
#Output: recipe_file path used
#Usage:  Setup recipe to use
sub setRecipeFile($) {
    my $recipe_file = shift;


    # if user provided recipe file
    if ( defined($recipe_file) ) {
#        my $tmp_recipe = $carun_cf->getOption((lc $recipe_file).'_recipe');
#
#        #First check if its a lookup in the conf file
#        if ( defined $tmp_recipe ) {
#            $recipe_file = $tmp_recipe;
#        }
#        #Else, instead of config lookup, a file is specified
#
    }
    else {
        $tf_object->bail("Recipe file not specified nor found in config file: $CARUN_CONFIG_FILE");
    }

    # Check if recipe file exists (from config or command-line)
    if ( !-e $recipe_file ) {
        $tf_object->bail("Recipe file does not exist: $recipe_file"
        );
    }

    $tf_object->logLocal( "Using recipe file as base: $recipe_file", 2 );

    return $recipe_file;
}

sub filterByStartEnd($$$) {
    my $recipe_file = shift;
    my $start = shift;
    my $end = shift;

    my $recipe_fh = new IO::File "<$recipe_file"
        or $tf_object->bail("Could not open recipe file: '$recipe_file'. Error code: $!");
    my @recipeLines = <$recipe_fh>;
    close $recipe_fh;

    ### ##########################################
    ### Set Start/End
    ### ##########################################
    #Number of original recipe lines
    my $size = @recipeLines;

    if ( defined $start ) {
        if ( $start =~ /^\d+$/){ #if digits provided
            $start--;    #Start is given by the user as 1-based.
        }
        else {    #if section name provided
            my $lineNumber = 0;
            foreach my $recipeLine (@recipeLines) {
              if ( $recipeLine =~ /^\s*#>\s*$start/) {
                $start = $lineNumber;
                last;
              }
              $lineNumber++;
            }
        }
    }
    else {
        $start = 0;
    }
    $tf_object->logLocal( "Using start: $start", 2 );

    if ( defined $end ) {
        if ( $end =~ /^\d+$/){ #if digits provided
            $end--;      #End is given by the user as 1-based
        }
        else {    #if section name provided
            my $lineNumber = 0;
            my $lineFound = 0;
            foreach my $recipeLine (@recipeLines) {
              if ( $lineFound == 1 ) {
                if ($recipeLine =~ /^\s*#>/) {
                    $end = $lineNumber - 1;
                    last;
                }
              }
              elsif ( $recipeLine =~ /#>\s*$end/) {
                  $lineFound = 1;
              }
              $lineNumber++;
            }
        }
    }
    else {
        $end = $size - 1;
    }
    $tf_object->logLocal( "Using end: $end", 2 );

    #The new set of lines to execute
    my @lines = @recipeLines[ $start .. $end ];

    #Open Temporary file for recipe
    my $timestamp = time();
    my $tmprecipe_file = "/tmp/.tarun.$timestamp";
    my $tmprecipe_fh = new IO::File ">$tmprecipe_file"
        or $tf_object->bail("Could not open temp recipe file: '$tmprecipe_file'. Error code: $!");
    foreach (@lines) {
        print $tmprecipe_fh $_;
    }
    close $tmprecipe_fh;
    $tf_object->logLocal( "Using temp recipe file as base: $tmprecipe_file", 2 );

    push @filesToClean, $tmprecipe_file;
    return $tmprecipe_file;
}

#Name:   passThroughSetup
#Input:  passThroughStr_ref
#        passThroughFile_ref
#Output: array of pass through commands
#Usage:  Setup Passthrough options
sub passThroughSetup ($) {
    my $passThroughStr_ref = shift;

    my @passThroughArr = ();
    #If passThroughStr is defined,
    #Arguments are space separated,
    #Each argument will put on a separate line
    #in the output script
    if ( defined ${$passThroughStr_ref} ) {
        @passThroughArr = split( / /, ${$passThroughStr_ref} );
    }
    return @passThroughArr;
}

sub generateSessionFile($) {
    my $caOptionsArray_ref = shift;

    my @taOptionsArray = @$caOptionsArray_ref;
    my $caOptionsStr = join('',@taOptionsArray);
    my $sessionFileName = "/tmp/options.sh";
    my $session_fh = new IO::File ">$sessionFileName"
      or $tf_object->bail("Could not open recipe file: '$sessionFileName'. Error Code: $!");
    print $session_fh "#TA Options:\n$caOptionsStr\n\n";
    close $session_fh;

    push @filesToClean, $sessionFileName;

    return $sessionFileName;
}

#Name:   taOptionsSetup
#Input:  erate
#        bubbleSmoothing
#        statMin
#        statMax
#        genomeSize
#Output: taOptionsArray - array of commands to run in the scripts
#Usage:  Checks for user provided TA options, and generates command lines for recipe
sub taOptionsSetup($) {
    my $Options           = shift;

    my @taOptionsArray = ();

    push @taOptionsArray, "OPTIONS=$Options\n";

    return @taOptionsArray;
}

sub validateErate($) {
    my $erate = shift;
    my $ERATE_LOW = $carun_cf->getOption('erate_low');
    my $ERATE_HIGH = $carun_cf->getOption('erate_high');

    $tf_object->bail("Invalid erate value '$erate'.  Valid values between $ERATE_LOW and $ERATE_HIGH")
        if ( $erate < $ERATE_LOW or $erate > $ERATE_HIGH);
}
use FindBin qw($Bin);
sub callArun($$$) {
    my $arun_options    = shift;
    my $sessionFileName = shift;
    my $recipeFileName  = shift;

    #my $arun = "$Bin/arun.pl";
    my $arun = "/usr/local/common/ARUN/arun";

    my $cmd = "cat  $sessionFileName $recipeFileName";
    $cmd .= "| $arun $arun_options" ;

    $tf_object->logLocal("Running '$cmd'",1);
    $tf_object->runCommand($cmd);
}

sub tigrFoundationOptions($$$$$$) {
    my $version     = shift;
    my $appendlog   = shift;
    my $logfile     = shift;
    my $debug       = shift;
    my $help        = shift;
    my $depend      = shift;

    $tf_object->printHelpInfoAndExit()
        if ( (defined $help) && ($help =~ /^(.*)$/) );

    $tf_object->printVersionInfoAndExit()
        if ( (defined $version) && ($version =~ /^(.*)$/) );

    $tf_object->printDependInfoAndExit()
        if ( (defined $depend) && ($depend =~ /^(.*)$/) );

    if ( (defined $appendlog) && ($appendlog =~ /^(.*)$/) ) {
       $appendlog = $1;
       $tf_object->logAppend($appendlog);
    }

    if ( (defined $logfile) && ($logfile =~ /^(.*)$/) )  {
       $logfile = $1;
       $tf_object->setLogFile($logfile);
    }

    if ( (defined $debug) && ($debug =~ /^(.*)$/) ) {
       $debug = $1;
       $tf_object->setDebugLevel($debug);
    }
}

sub clean() {
    foreach my $file (@filesToClean) {
        unlink $file or $tf_object->bail("Unable to remove: $file");
    }
}

#Name:   SIGHANDLER
#Input:  none
#Output: none
#Usage:  this functions acts as the signal handler for the external signals
#        that caexec might receive. It indicates that an external signal has
#        been received by changing the status field in the Request table to
#        "C" (cancel).
sub SIGHANDLER {
    $tf_object->logLocal( "tarun received a cancel request", 1 );
    print STDOUT "tarun received a cancel request\n";
    $cancel_track = 1;
    exit(2);
}

# =============================== MAIN ======================================
#
MAIN:
{

 #################################### BEGIN: borrowed from carun ################################
    my $recipe_file     = undef;    # Recipe to run
    my $passThroughStr  = undef;    # passthrough commands (via command line)

    #TIGR Foundation variables
    my $version         = undef;    # Version option
    my $help            = undef;
    my $depend          = undef;
    my $appendlog       = undef;
    my $logfile         = undef;
    my $debug           = undef;
    my $alias           = undef;

    # ========================== Program Setup ==============================
    # Prepare logs
    $tf_object->addDependInfo(@MY_DEPENDS);
    $tf_object->setHelpInfo($HELPTEXT);
    $tf_object->setVersionInfo($MY_VERSION);

    Getopt::Long::Configure('pass_through');
    Getopt::Long::Configure('no_ignore_case');
    GetOptions(
        'recipe=s',       \$recipe_file,
        'p=s',       \$passThroughStr,
        'version|V', \$version,         'appendlog=i', \$appendlog,
        'logfile=s', \$logfile,         'debug=i', \$debug,
        'help|h',    \$help,            'depend', \$depend,
        'alias=s',   \$alias
    ) or $tf_object->bail("Options could not be read. See tarun -h.");


    tigrFoundationOptions($version, $appendlog, $logfile, $debug, $help, $depend);

 #################################### END: borrowed from carun ##################################
  $debug            = 0;
  my $destination = undef;
  my $local = 0;
  my $bigjob = 0;

  # Command line parms control TIGR_Assembler
  my $UseConstraints;      undef $UseConstraints;
  my $UseCloneConstraints; undef $UseCloneConstraints;
  my $UseMultiMatchRecs;   undef $UseMultiMatchRecs;
  my $SeqFile;             undef $SeqFile;
  my $ContigFile;          undef $ContigFile;
  my $QualFile;            undef $QualFile;
  my $MaxEnd;              undef $MaxEnd;
  my $MaxMatchLength;      undef $MaxMatchLength;
  my $MinLen;              undef $MinLen;
  my $MinPct;              undef $MinPct;
  my $MaxErr32;            undef $MaxErr32;
  my $UseSingletons;       undef $UseSingletons;
  my $StartPct;            undef $StartPct;
  my $EndPct;              undef $EndPct;
  my $PctStep;             undef $PctStep;
  my $SafeFileSize;        undef $SafeFileSize;
  # Command line parms for test
  my $tigrasm          = "";


  # ========================== Program Setup ==============================
  $tf_object->setDebugLevel(0);
  my $basedir = $tf_object->getProgramInfo("env_path");

  # now we handle the input options
  my $result = undef;
  if ($MY_APP eq "tarun")
  {
    $result = $tf_object->TIGR_GetOptions
              (
                 # Control tarun
                 'local',    \$local,
                 'dir:s',    \$destination,
	         	 'bigjob',   \$bigjob,
	         	 'recipe=s', \$recipe_file,

                 # Control the TIGR_Assembler
                 'constraints!',   \$UseConstraints,
                 'C:s',            \$ContigFile,
                 'e:s',            \$MaxEnd,
                 'j!',             \$UseCloneConstraints,
                 'g:i',            \$MaxErr32,
                 'l:i',            \$MinLen,
                 'M!',             \$UseMultiMatchRecs,
                 'p:f',            \$MinPct,
                 'q:s',            \$QualFile,
                 's!',             \$UseSingletons,
                 'Y:i',            \$SafeFileSize,
                 # Undocumented TIGR_Assembler options (pass through)
                 'a:s', \&appOptions, # alignment_directory
                 'A:s', \&appOptions, # ace_output_file
                 'b:i', \&appOptions, # buffer size
                 'c:s', \&appOptions, # coverage_filename
                 'd',   \&appFlags,   # -
                 'D:s', \&appOptions, # phd_dir
                 'E:i', \&appOptions, # phd_dir
                 'f:s', \&appOptions, # fasta_file
                 'F:s', \&appOptions, # repeat_file
                 'G:s', \&appOptions, # grouper file
                 'J:s', \&appOptions, # dump_file
                 'L',   \&appFlags,   # -
                 'n:s', \&appOptions, # asm_prefix
                 'N',   \&appFlags,   # -
                 'o:s', \&appOptions, # score file
                 'P',   \&appFlags,   # -
                 'r:i', \&appOptions, # resort number
                 'R',   \&appFlags,   # -
                 'S:i', \&appOptions, # max_span_len
                 't',   \&appFlags,   # -
                 'T',   \&appFlags,   # -
                 'u',   \&appFlags,   # -
                 'w:i', \&appFlags,   # number of threads
                 'x',   \&appFlags,   # -
                 'X',   \&appFlags,   # -
                 'y:i', \&appOptions, # repeat_num_cutoff
                 'z:i', \&appOptions, # num_conflicts
                 'Z',   \&appFlags,   # -
                 # Test Controls
                 'tigrasm:s',      \$tigrasm,
              );
  }
  else
  {
    $result = $tf_object->TIGR_GetOptions
              (
                 # Control tarun_beta
                 'local',     \$local,
                 'dir:s',     \$destination,
	         'bigjob',    \$bigjob,
	         'recipe=s',       \$recipe_file,

                 # Control the TIGR_Assembler
                 'constraints!',   \$UseConstraints,
                 'C:s',            \$ContigFile,
                 'e:s',            \$MaxEnd,
                 'j!',             \$UseCloneConstraints,
                 'k:i',            \$MaxMatchLength,
                 'g:i',            \$MaxErr32,
                 'l:i',            \$MinLen,
                 'M!',             \$UseMultiMatchRecs,
                 'p:f',            \$MinPct,
                 'q:s',            \$QualFile,
                 's!',             \$UseSingletons,
                 'U:f',            \$StartPct,
                 'B:f',            \$EndPct,
                 'i:f',            \$PctStep,
                 'Y:i',            \$SafeFileSize,
                 # Undocumented TIGR_Assembler options (pass through)
                 'a:s', \&appOptions, # alignment_directory
                 'A:s', \&appOptions, # ace_output_file
                 'b:i', \&appOptions, # buffer size
                 'c:s', \&appOptions, # coverage_filename
                 'd',   \&appFlags,   # -
                 'D:s', \&appOptions, # phd_dir
                 'E:i', \&appOptions, # phd_dir
                 'f:s', \&appOptions, # fasta_file
                 'F:s', \&appOptions, # repeat_file
                 'G:s', \&appOptions, # grouper file
                 'J:s', \&appOptions, # dump_file
                 'L',   \&appFlags,   # -
                 'n:s', \&appOptions, # asm_prefix
                 'N',   \&appFlags,   # -
                 'o:s', \&appOptions, # score file
                 'P',   \&appFlags,   # -
                 'r:i', \&appOptions, # resort number
                 'R',   \&appFlags,   # -
                 'S:i', \&appOptions, # max_span_len
                 't',   \&appFlags,   # -
                 'T',   \&appFlags,   # -
                 'u',   \&appFlags,   # -
                 'w:i', \&appFlags,   # number of threads
                 'x',   \&appFlags,   # -
                 'X',   \&appFlags,   # -
                 'y:i', \&appOptions, # repeat_num_cutoff
                 'z:i', \&appOptions, # num_conflicts
                 'Z',   \&appFlags,   # -
                 # Test Controls
                 'tigrasm:s',      \$tigrasm,
              );
  }
  $tf_object->bail("Command line parsing failed") if ($result == 0);

  # Test settings
  # Use instead the specified TIGR Assembler (hidden test option)
  if ($tigrasm ne "")
  {
    $TIGRASM = $tigrasm;
    $LOCAL_TIGRASM = "$tigrasm";  # rebind this value
  }
  # Pass-through some TIGR_Assembler options.  All are single letter,
  # may contain additional alphanumeric field, including pathname


  # Resolve relative paths and check for accessibility of input files.
  # There must be at least a sequence file on input.
  $SeqFile = $ARGV[-1];   # Should be leftover after GetOpts processing
  $tf_object->bail("Missing .seq file ...") if (! defined $SeqFile);

  my ($name, $path, $suffix) = fileparse($SeqFile, ".seq");
  my $outfile = $name;

  # Create canonical paths for condor
  $ContigFile = "$basedir/$ContigFile" if (defined $ContigFile);
  $QualFile   = "$basedir/$QualFile"   if (defined $QualFile);
  $SeqFile   = "$basedir/$SeqFile"   if (defined $SeqFile);

  # In 2.5 beta and later a safe resource allocation feature ensures that
  # the correct amount of space is computed in the beginning of the assembly
  # so it doesn't run out of space late in the assembly.
  if ($MY_SAFEOPS && ($SeqFile ne ""))
  {
    # Compute it if it's not given at the command line
    if (! defined $SafeFileSize)
    {
      # Get the size of the .seq file
      my $file_size = (stat $SeqFile)[7];
      $tf_object->bail("Cannot access sequence file \'$SeqFile\'") if (! defined $file_size);

      # Find the percent value used.
      my $percent_value = (defined $MinPct)? $MinPct : $DEFAULT_MIN_PCT;
      my $percent_expected_gaps = (101 - $percent_value) / 100;
      my $gap_overhead = int ($file_size * $percent_expected_gaps);
      $SafeFileSize = $file_size + $gap_overhead;
    }
  }

  # Prepare the list of options for TIGR_Assembler from internal parms.  Append
  # the undocumented pass-through options to the end of this list in case some of
  # them override the existing defaults.  Translation of use-case oriented options
  # into TA options (ie -constraints) is done here.
  my $Options = "";
  $Options .= "-j "      if (($UseConstraints && (! defined $UseCloneConstraints)) ||
                             (defined $UseCloneConstraints && $UseCloneConstraints)
                            );
  $Options .= "-M "      if (($UseConstraints && (! defined $UseMultiMatchRecs))   ||
                             (defined $UseMultiMatchRecs   && $UseMultiMatchRecs)
                            );
  $Options .= "-C $ContigFile "      if (defined $ContigFile);
  $Options .= "-e $MaxEnd "          if (defined $MaxEnd);
  $Options .= "-k $MaxMatchLength "  if (defined $MaxMatchLength);
  $Options .= "-g $MaxErr32 "        if (defined $MaxErr32);
  $Options .= "-l $MinLen "          if (defined $MinLen);
  $Options .= "-p $MinPct "          if (defined $MinPct);
  $Options .= "-q $QualFile "        if (defined $QualFile);
  $Options .= "-s "                  if (defined $UseSingletons && $UseSingletons);
  $Options .= "-U $StartPct "        if (defined $StartPct);
  $Options .= "-B $EndPct "          if (defined $EndPct);
  $Options .= "-i $PctStep "         if (defined $PctStep);
  $Options .= "-Y $SafeFileSize "    if (defined $SafeFileSize);
  $Options .= $APP_OPTIONS;

  # Prepare the output directory.  By default it is $PWD/asm_<timestamp>
  # In the future output objects will live in a directory stucture that will be
  # associated with a revision system.  This needs to be 777 on condor.
  #
  my ($sec, $min, $hour, $day, $mon, $year, @rest) = localtime(time);
  my $mydate = sprintf($MY_DATE_FORMAT, $year+1900, $mon+1, $day);
  my $outdir = (defined $destination)? $destination : "$basedir/asm_$mydate";
  # Find a unique directory name so that concurrent execution is possible.
  my $chr = 'a';
  while (-e $outdir)
  {
    $outdir = (defined $destination)? $destination : "$basedir/asm_$mydate";
    $outdir .= $chr;
    $chr++;
  }
  umask(0000);
  mkdir($outdir, 0777);
  if ( !(-d $outdir)  ||  !(-x $outdir)  ||  !(-w $outdir) )
  {
    $tf_object->bail("Could not access output directory $outdir.");
  }


  # Prepare the filenames
  my $scratch_file = "$outdir/$outfile.scratch";
  my $asm_file     = "$outdir/$outfile.tasm";
  my $align_file   = "$outdir/$outfile.align";
  my $fasta_file   = "$outdir/$outfile.fasta";
  my $error_file   = "$outdir/$outfile.error";
  my $log_file     = "$outdir/$outfile.log";
#  my $store_file   = ($MY_DISKCACHE) ?  "$outdir/outfile.store" : "";

  if ($local)
  {
    # Check the environment
    $tf_object->bail("Could not run local assembly ($LOCAL_TIGRASM not executable).", 1)
      if (! -e $LOCAL_TIGRASM  ||  ! -x $LOCAL_TIGRASM);

    # Prepare the options and launch it.
    my $cmd = "$LOCAL_TIGRASM";
    $cmd .= " $Options -n $outfile -a $align_file -f $fasta_file $scratch_file " .
            " < $SeqFile > $asm_file 2>$error_file";
    $tf_object->logLocal("Running $cmd ...", 0);
    my $result = $tf_object->runCommand($cmd);
    if ($result != 0)
    {
      # Remove the scratch file if it exists
      if (-e $scratch_file)
      {
        unlink $scratch_file;  # Don't report failure
      }
      $tf_object->bail("Command failed: $cmd ($result)");
    }
  }
  else
  {
    # Check the environment
    my $OPTIONS = "\$OPTIONS=".$Options;
    my @f = stat $SeqFile;

    my $fastoption = "-l fast";
    my $passthru = ($bigjob) ? "" : "--passthrough \"$fastoption\" " ;
#    my $run_command = "qsub   -M $user -m e \"$assembler_cmd\" "."  $passthru ";
    # Set recipe file
    $recipe_file = setRecipeFile($recipe_file);

    # Set option values
    my @taOptionsArray = taOptionsSetup("\"".$Options."\"");
#   my @passThroughCommands = passThroughSetup (\$passthru);
    my $sessionFileName = generateSessionFile(\@taOptionsArray);
    my $invocation = "$0 ".$tf_object->getProgramInfo("invocation");

    my $arunOptions = join(' ', @ARGV);
    $arunOptions .= " -invocation \"$invocation\"";
    $arunOptions .= " -alias \"$alias\"" if (defined $alias);

    # Call Arun
    callArun($arunOptions,$sessionFileName,$recipe_file);
#    callArun($arunOptions,$recipe_file);
    # Clean up
    clean();
  }
  exit 0;
}
