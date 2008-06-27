#!/usr/local/bin/perl -w
# $Id: carun.pl,v 1.4 2008-06-27 06:29:19 brianwalenz Exp $
#
# Celera Assembler front-end script calls runs assembly pipline
# recipes on the grid or on the local host.
#
# Written by Marwan Oweis September 2006.
#

# Program configuration
my @MY_DEPENDS = ( "TIGR::Foundation", "TIGR::ConfigFile" );
my $MY_VERSION = " 1.5 (Build " . (qw/$Revision: 1.4 $/)[1] . ")";
my $HELPTEXT =
  qq~Request a whole-genome shotgun assembly using Celera Assembler.

 usage: carun [options] <frg>
	    carun -cancel <id>

 inputs:
  <frg>       Sequence file generated e.g. by pullfrag.

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
  -R <recipe>       Specify recipe to execute Celera Assembler
                    [no]OBT - run with OBT
                    noVec   - run w/OBT but w/o Vector
                    hybrid  - run w/OBT using hybrid reads
  -before           Print the recipe and exit
  -after            Print the generated shell script and exit

 CA options:
  -e <percent>      Unitigger assumed percent error [0.00,0.06] (default: 0.015)
  -g <len>          Assign genome length for unitigger
  -j <lbound>       Set scaffolder A-stat low bound (default: 1)
  -k <hbound>       Set scaffolder A-stat high bound (default: 5)
  -[no]ubs          Enable unitigger bubble smoothing (default: enabled)

 advanced options:
  -config <file>    Specify config file to be used
  -cont <req_id>    Use the request id given for continuation.
  -start <d>        Start at line <d> in recipe file. (1-based, Inclusive)
                    Alternatively, section names can be used.  Section names
                    begin with '#>' in a recipe file.
  -end  <d>         End at line <d> in recipe file. (1-based, Inclusive)
                    Alternatively, section names can be used.  Section names
                    begin with '#>' in a recipe file.
  -p <args>         Passthrough arguments inserted at beginning of recipe,
                    each argument on a separate line (space-delimited).
  -P <file>         Same as -p, but retrieved from file.  If both specified,
                    -p arguments will be inserted after file-retrieved content
                    (command-line has precendence over file-retrieved content).
  -test             Run in debug mode (same as -D test -nonotify)

 Genomic Assembly Pipeline:
 https://intranet.jcvi.org/cms/SE/GAP
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

# ================================== Constants  ===============================
# default config file for carun
my $DEFAULT_INSTALL_DIR = "/usr/local/common/ARUN";
my $CARUN_CONFIG_FILE = "CARUN.conf";
my $ARUN_CONFIG_FILE = "ARUN.conf";
my $install_dir = $ENV{'ARUN_INSTALL_DIR'};

# The carun Config file object
my $carun_cf = undef;
my $arun_cf  = undef;
my $timestamp = time();

# when Cntl-C is used to cancel carun
my $cancel_track = 0;

# running at JCVI
my $jcvi = undef;

# Module for housekeeping functions
my $tf_object = new TIGR::Foundation;

# Stack of files to clean at the end of the run
my @filesToClean = ();

#Name:   initializeConfig
#Input:  config file name
#Usage:  Retrive config file values
sub initializeConfig($) {
    my $configFile = shift;

    $install_dir = $DEFAULT_INSTALL_DIR if ( !defined $install_dir );
    $configFile = $install_dir . '/' . $CARUN_CONFIG_FILE if ( !defined $configFile );
    my $arunConfigFile = $install_dir . '/' . $ARUN_CONFIG_FILE;

    # The carun Config file object
    $carun_cf = new TIGR::ConfigFile($configFile)
      or $tf_object->bail("Could not initialize the config file: '$configFile'");
    $tf_object->logLocal( "Using carun config file: $configFile", 2 );

    # The arun Config file object
    $arun_cf = new TIGR::ConfigFile($arunConfigFile)
      or $tf_object->bail("Could not initialize the config file: '$arunConfigFile'");
    $tf_object->logLocal( "Using arun config file: $arunConfigFile", 2 );
}


#Name:   setRecipeFile
#Input:  recipe_file per user
#Output: recipe_file path used
#Usage:  Setup recipe to use
sub setRecipeFile($$$) {
    my $recipe_file = shift;
    my $start = shift;
    my $end = shift;

    #default_recipe option has the name of the default recipe option
    my $default_recipe_file =
      $carun_cf->getOption( $carun_cf->getOption('default_recipe') );

    # if user provided recipe file
    if ( defined($recipe_file) ) {
        my $tmp_recipe = $carun_cf->getOption((lc $recipe_file).'_recipe');

        #First check if its a lookup in the conf file
        if ( defined $tmp_recipe ) {
            $recipe_file = $install_dir . '/'  . $tmp_recipe;
        }
        #Else, instead of config lookup, a file is specified

    }
    elsif ( defined($default_recipe_file) ) {
        $recipe_file = $install_dir . '/'  . $default_recipe_file;
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

    $recipe_file = filterByStartEnd($recipe_file, $start, $end)
        if ((defined $start) or (defined $end));

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
    my $tmprecipe_file = "/tmp/carun.$timestamp.recipe";
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
sub passThroughSetup ($$) {
    my $passThroughStr_ref = shift;
    my $passThroughFile_ref = shift,

    my @passThroughArr = ();
    #If passThroughStr is defined,
    #Arguments are space separated,
    #Each argument will put on a separate line
    #in the output script
    if ( defined ${$passThroughStr_ref} ) {
        @passThroughArr = split( / /, ${$passThroughStr_ref} );
    }

    if ( defined ${$passThroughFile_ref} ) {
        $tf_object->bail("Passthrough file: ${$passThroughFile_ref} is not found.")
          if ( !-e ${$passThroughFile_ref} );
        my $passThru_fh = new IO::File "<${$passThroughFile_ref}"
          or $tf_object->bail("Could not open pass through file: '${$passThroughFile_ref}'. Error code: $!");
        my (@lines) = <$passThru_fh>;
        close $passThru_fh
            or $tf_object->bail("Error closing '${$passThroughFile_ref}'. Error Code: $!");

        #Add file arguments at the beginning of the passThroughArr
        #Command-line arguments have precedence over file arguments,
        #hence they will come after as to overwrite (if applicable)
        #the file arguments
        unshift( @passThroughArr, @lines );
    }

    return @passThroughArr;
}

sub setPostRecipeFile () {
    $tf_object->logLocal( "Inside setPostRecipeFile", 1 );

    my $pr_str = q~
#> Upload Metrics
cd $WORKDIR
/local/packages/aserver/ametrics.plx -debug 9 $REQUEST_ID $PREFIX.qc.metrics

#> Cleanup
rm $WORKDIR/$PREFIX.cgi $WORKDIR/$PREFIX.frg
~;

    my $prFileName = "/tmp/carun.$timestamp.pr.sh";
    my $pr_fh = new IO::File ">$prFileName"
      or $tf_object->bail("Could not create post recipe file: '$prFileName'. Error Code: $!");
    print $pr_fh $pr_str;
    close $pr_fh;

    push @filesToClean, $prFileName;

    return $prFileName;

}
sub generateSessionFile($$) {
    my $caOptionsArray_ref = shift;
    my $passThroughCommands_ref = shift;

    my @caOptionsArray = @$caOptionsArray_ref;
    my @passThroughArray = @$passThroughCommands_ref;
    my $caOptionsStr = join('',@caOptionsArray);
    my $passThroughStr = join('',@passThroughArray);
    my $sessionFileName = "/tmp/carun.$timestamp.options.sh";
    my $session_fh = new IO::File ">$sessionFileName"
      or $tf_object->bail("Could not open recipe file: '$sessionFileName'. Error Code: $!");
    print $session_fh "#CA Options:\n$caOptionsStr\n\n";
    print $session_fh "#Pass Through:\n$passThroughStr\n\n";
    close $session_fh;

    push @filesToClean, $sessionFileName;

    return $sessionFileName;
}

#Name:   caOptionsSetup
#Input:  erate
#        bubbleSmoothing
#        statMin
#        statMax
#        genomeSize
#Output: caOptionsArray - array of commands to run in the scripts
#Usage:  Checks for user provided CA options, and generates command lines for recipe
sub caOptionsSetup($$$$$) {
    my $erate           = shift;
    my $bubbleSmoothing = shift;
    my $statMin         = shift;
    my $statMax         = shift;
    my $genomeSize      = shift;

    my @caOptionsArray = ();

    # setup recipe specific parameters
    $erate = $carun_cf->getOption('erate_default')
        if ( !defined $erate );

    validateErate($erate);
    $tf_object->logLocal( "Using erate=$erate%", 2 );

    push @caOptionsArray, "ERATE=$erate\n";
    push @caOptionsArray, "echo ERATE=\$ERATE\n";

    if ( defined $bubbleSmoothing && $bubbleSmoothing eq 0 ) {
        $bubbleSmoothing = $carun_cf->getOption('bubble_disable');
        $tf_object->logLocal( "Using Bubble Smoothing=disabled", 2 );
    }
    else {
        $bubbleSmoothing = $carun_cf->getOption('bubble_enable');
        $tf_object->logLocal( "Using Bubble Smoothing=enabled", 2 );
    }
    push @caOptionsArray, "BUBBLE=$bubbleSmoothing\n";
    push @caOptionsArray, "echo BUBBLE=\$BUBBLE\n";

    $statMin = $carun_cf->getOption('statlow_default')
        if ( !defined $statMin );
    $tf_object->logLocal( "Using Astat Low Bound=$statMin", 2 );
    push @caOptionsArray, "ASTATLOW=$statMin\n";
    push @caOptionsArray, "echo ASTATLOW=\$ASTATLOW\n";

    $statMax = $carun_cf->getOption('stathigh_default')
        if ( !defined $statMax );
    $tf_object->logLocal( "Using Astat High Bound=$statMax", 2 );
    push @caOptionsArray, "ASTATHIGH=$statMax\n";
    push @caOptionsArray, "echo ASTATHIGH=\$ASTATHIGH\n";

    if ( defined $genomeSize ) {
        push @caOptionsArray, "GENOMELENGTH_L=\"-l $genomeSize\"\n";
        push @caOptionsArray, "echo GENOMELENGTH_L=\$GENOMELENGTH_L\n";
        push @caOptionsArray, "GENOMELENGTH_G=\"-g $genomeSize\"\n";
        push @caOptionsArray, "echo GENOMELENGTH_G=\$GENOMELENGTH_G\n";
    }

    $jcvi = $arun_cf->getOption('jcvi');
    if( $jcvi == 1) {
        push @caOptionsArray, "EUIDSERVICE=\"-u -E \$EUIDList -B \$EUIDBlockSize -n \$EUIDNamespace\"\n";
        push @caOptionsArray, " echo EUIDSERVICE=\$EUIDSERVICE\n";
    }

    return @caOptionsArray;
}

sub validateErate($) {
    my $erate = shift;
    my $ERATE_LOW = $carun_cf->getOption('erate_low');
    my $ERATE_HIGH = $carun_cf->getOption('erate_high');

    $tf_object->bail("Invalid erate value '$erate'.  Valid values between $ERATE_LOW and $ERATE_HIGH")
        if ( $erate < $ERATE_LOW or $erate > $ERATE_HIGH);
}

sub callArun($$$$$) {
    my $arun_options    = shift;
    my $sessionFileName = shift;
    my $recipeFileName  = shift;
    my $postRecipeFileName = shift;
    my $before          = shift;

    my $carunSetup = $install_dir . '/'  . $carun_cf->getOption('carun_setup');
    my $arun = $install_dir . '/'  . $carun_cf->getOption('arun_exec');

    my $cmd = "cat $carunSetup $sessionFileName $recipeFileName $postRecipeFileName";
    $cmd .= "| $arun $arun_options" if (!$before);

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
    } else {
       $tf_object->setDebugLevel($carun_cf->getOption('debug_default'));
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
    $tf_object->logLocal( "carun received a cancel request", 1 );
    print STDOUT "carun received a cancel request\n";
    $cancel_track = 1;
    exit(2);
}

# =============================== MAIN ======================================
#
MAIN:
{
    my $recipe_file     = undef;    # Recipe to run
    my $post_recipe_file= undef;    # To run after the recipe
    my $before          = undef;    # Output recipe file
    my $erate           = undef;    # CA option
    my $genomeSize      = undef;    # CA option
    my $bubbleSmoothing = undef;    # CA option
    my $statMin         = undef;    # CA option
    my $statMax         = undef;    # CA option
    my $configFile      = undef;    # optional user specified config file
    my $passThroughStr  = undef;    # passthrough commands (via command line)
    my $passThroughFile = undef;    # passthrough commands (via file)
    my $start           = undef;    # user specified line number in recipe to start
    my $end             = undef;    # user specified line number in recipe to end

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
        'R=s',       \$recipe_file,     'before',    \$before,
        'e=f',       \$erate,           'g=i',       \$genomeSize,
        'ubs!',      \$bubbleSmoothing, 'j=i',       \$statMin,
        'k=i',       \$statMax,         'config=s',  \$configFile,
        'p=s',       \$passThroughStr,  'P=s',       \$passThroughFile,
        'start=s',   \$start,           'end=s',     \$end,
        'version|V', \$version,         'appendlog=i', \$appendlog,
        'logfile=s', \$logfile,         'debug=i', \$debug,
        'help|h',    \$help,            'depend', \$depend,
        'alias=s',   \$alias
    ) or $tf_object->bail("Options could not be read. See carun -h.");

    # Initialize parameters based on config file
    initializeConfig($configFile);

    tigrFoundationOptions($version, $appendlog, $logfile, $debug, $help, $depend);

    # Set recipe file
    $recipe_file = setRecipeFile($recipe_file,$start,$end);

    # Set option values
    my @caOptionsArray = caOptionsSetup( $erate, $bubbleSmoothing, $statMin, $statMax, $genomeSize);
    my @passThroughCommands = passThroughSetup (\$passThroughStr, \$passThroughFile);
    my $sessionFileName = generateSessionFile(\@caOptionsArray,\@passThroughCommands);
    my $invocation = "$0 ".$tf_object->getProgramInfo("invocation");

    # Set post recipe file
    $post_recipe_file = '';
    $post_recipe_file = setPostRecipeFile() if ( $jcvi == 1 );

    my $arunOptions = join(' ', @ARGV);
    $arunOptions .= " -invocation \"$invocation\"";
    $arunOptions .= " -alias \"$alias\"" if (defined $alias);

    # Call Arun
    callArun($arunOptions,$sessionFileName,$recipe_file,$post_recipe_file,$before);

    # Clean up
    clean();
}
