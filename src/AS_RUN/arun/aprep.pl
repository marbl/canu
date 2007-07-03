#!/usr/local/bin/perl -w
# $Id: aprep.pl,v 1.1 2007-07-03 19:44:35 moweis Exp $
#
# Celera Assembler front-end script calls runs assembly pipline 
# recipes on the grid or on the local host.
#
# Written by Marwan Oweis September 2006.
#

# Program configuration
my @MY_DEPENDS = ( "TIGR::Foundation", "TIGR::ConfigFile" );
my $MY_VERSION = " 1.2 (Build " . (qw/$Revision: 1.1 $/)[1] . ")";
my $HELPTEXT =
  qq~Perform preload processes on a Celera Assembly.

 usage: aprep [options] <prefix>
	    aprep -cancel <id>
 
 inputs: 
  <prefix>       Assembly prefix

 general options:
  -alias <a>        Identify the request as <a> on Assembly Server Console
  -cancel <id> 	    Cancel assembly request <id>
  -noCopy           Do not copy. (Default: copy to 'aprepCopy')
  -copyDir <dir>    Specify dir to copy to (instead of 'aprepCopy')
  -offgrid, -local  Assemble on this host (default: on grid)
  -D <name>         Specify project name (displayed on AserverConsole).
                    Use -D test for debugging. (default: frg prefix) 
  -[no]notify    	Send email upon completion/error and create BITS case 
                    for errors. (default: notify)

 recipe options:
  -R <recipe>       Specify recipe to execute. 
  -before           Print the recipe and exit
  -after            Print the generated shell script and exit
 
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
      
 Genomic Assembly Pipeline: 
 https://intranet.jcvi.org/cms/SE/GAP
 Tracking assembly requests:  
 http://aserver.tigr.org:8080/AserverConsole/
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
# default config file for aprep
my $DEFAULT_INSTALL_DIR = "/usr/local/common/CARUN";
my $PREP_CONFIG_FILE = "APREP.conf";
my $install_dir = $ENV{'ARUN_INSTALL_DIR'};    

# The aprep Config file object
my $prep_cf = undef;

# when Cntl-C is used to cancel aprep
my $cancel_track = 0;

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
    $configFile = $install_dir . '/' . $PREP_CONFIG_FILE if ( !defined $configFile );

    # The aprep Config file object
    $prep_cf = new TIGR::ConfigFile($configFile)
      or $tf_object->bail("Could not initialize the config file: '$configFile'");
    $tf_object->logLocal( "Using config file: $configFile", 2 );    
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
      $prep_cf->getOption( $prep_cf->getOption('default_recipe') );

    # if user provided recipe file
    if ( defined($recipe_file) ) {
        my $tmp_recipe = $prep_cf->getOption((lc $recipe_file).'_recipe');

        #First check if its a lookup in the conf file
        if ( defined $tmp_recipe ) {
            $recipe_file = $install_dir . '/' . $tmp_recipe;
        }
        #Else, instead of config lookup, a file is specified
        
    }
    elsif ( defined($default_recipe_file) ) {
        $recipe_file = $install_dir . '/' . $default_recipe_file;
    }
    else {
        $tf_object->bail("Recipe file not specified nor found in config file: $PREP_CONFIG_FILE");
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
    my $timestamp = time();
    my $tmprecipe_file = "/tmp/.aprep.$timestamp";
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

sub generateSessionFile($) {    
    my $passThroughCommands_ref = shift;
    
    my @passThroughArray = @$passThroughCommands_ref;
    my $passThroughStr = join('',@passThroughArray);     
    my $sessionFileName = "/tmp/options.sh";
    my $session_fh = new IO::File ">$sessionFileName"    
      or $tf_object->bail("Could not open recipe file: '$sessionFileName'. Error Code: $!");    
    print $session_fh "#Pass Through:\n$passThroughStr\n\n";    
    close $session_fh;    
    push @filesToClean, $sessionFileName;    

    return $sessionFileName;         
}

sub callArun($$$$) {
    my $arun_options    = shift;
    my $sessionFileName = shift;
    my $recipeFileName  = shift;
    my $before          = shift;

    my $prepSetup = $install_dir . '/' . $prep_cf->getOption('prep_setup');
    my $arun = $install_dir . '/' . $prep_cf->getOption('arun_exec');

    my $cmd = "cat $prepSetup $sessionFileName $recipeFileName";
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
       $tf_object->setDebugLevel($prep_cf->getOption('debug_default'));
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
    $tf_object->logLocal( "aprep received a cancel request", 1 );
    print STDOUT "aprep received a cancel request\n";
    $cancel_track = 1;
    exit(2);
}

# =============================== MAIN ======================================
#
MAIN:
{
    my $recipe_file     = undef;    # Recipe to run
    my $before          = undef;    # Output recipe file
    my $configFile      = undef;    # optional user specified config file
    my $passThroughStr  = undef;    # passthrough commands (via command line)    
    my $passThroughFile = undef;    # passthrough commands (via file)    
    my $start           = undef;    # user specified line number in recipe to start
    my $end             = undef;    # user specified line number in recipe to end
    my $alias           = undef;
    my $copyDir         = undef;
    my $noCopy          = undef;

    #TIGR Foundation variables
    my $version         = undef;    # Version option
    my $help            = undef;
    my $depend          = undef;
    my $appendlog       = undef;
    my $logfile         = undef;
    my $debug           = undef;
    
    # ========================== Program Setup ==============================
    # Prepare logs
    $tf_object->addDependInfo(@MY_DEPENDS);
    $tf_object->setHelpInfo($HELPTEXT);
    $tf_object->setVersionInfo($MY_VERSION);

    Getopt::Long::Configure('pass_through');
    Getopt::Long::Configure('no_ignore_case');
    GetOptions(
        'R=s',       \$recipe_file,     'before',    \$before,
        'config=s',  \$configFile,      'alias=s',   \$alias,                        
        'p=s',       \$passThroughStr,  'P=s',       \$passThroughFile,
        'start=s',   \$start,           'end=s',     \$end,        
        'version|V', \$version,         'appendlog=i', \$appendlog,
        'logfile=s', \$logfile,         'debug=i', \$debug,
        'help|h', \$help,               'depend', \$depend,
        'copyDir=s', \$copyDir,         'noCopy', \$noCopy
    ) or $tf_object->bail("Options could not be read. See aprep -h.");

    # Initialize parameters based on config file
    initializeConfig($configFile);

    tigrFoundationOptions($version, $appendlog, $logfile, $debug, $help, $depend);
    
    bail("Please choose either copy <dir> or noCopy.") if ( defined $copyDir and defined $noCopy );
    
    # Set recipe file
    $recipe_file = setRecipeFile($recipe_file,$start,$end);

    # Set option values
    my @passThroughCommands = passThroughSetup (\$passThroughStr, \$passThroughFile);
    my $sessionFileName = generateSessionFile(\@passThroughCommands);

    # Call Arun    
    my $invocation = "$0 ".$tf_object->getProgramInfo("invocation");
    
    my $arunOptions = join(' ', @ARGV);
    $arunOptions .= " -invocation \"$invocation\"";
    $arunOptions .= " -alias \"$alias\"" if (defined $alias);
    
    my $aget_file = $install_dir . '/' . $prep_cf->getOption('aget_file');
    
    if ( !defined $noCopy ) {
        $arunOptions .= " -custCopy $aget_file";        
        $arunOptions .= ",$copyDir" if ( defined $copyDir and $copyDir ne '');
    }

    callArun($arunOptions,$sessionFileName,$recipe_file,$before);

    # Clean up
    clean();
}
