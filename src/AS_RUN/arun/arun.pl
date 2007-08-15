#!/usr/local/bin/perl
# $Id: arun.pl,v 1.2 2007-08-15 22:30:50 moweis Exp $
#
# Given input from STDIN, generates and runs shell script
# on the grid or on the local machine with AserverConsole support.
#
# Written by Marwan Oweis September 2006.
#

# Program configuration
my @MY_DEPENDS = ( "TIGR::Foundation", "TIGR::ConfigFile" );
my $MY_VERSION = " 1.7 (Build " . (qw/$Revision: 1.2 $/)[1] . ")";
my $HELPTEXT =
  qq~
Given input from STDIN, generates and runs shell script
on the grid or on the local machine with Aserver Console support.

 usage: arun [options] <prefix> [copy_options]
	    arun -cancel <id>
 
 copy options (only one of the following options allowed) (default: medCopy):
  -maxCopy [<dir>]  Copy max output to <dir> (default: current directory)
  -medCopy [<dir>]  Copy med output to <dir> (default: current directory) 
  -minCopy [<dir>]  Copy min output to <dir> (default: current directory)
  -custCopy <f>[,<dir>] Customized copy.  Two parameters given (comma-separated):
                    -File name: lists the file extentions to copy
                    -Output dir name: name of local directory to copy to 
                    (default: current directory)
  -noCopy           Do not make local copy.
  
 options:
  -alias <a>        Identify the request as <a> on Aserver Console
  -cancel <id> 	    Cancel assembly request <id>
  -config <file>    Specify config file to be used
  -fast             To bypass grid queue.  Recommended usage for small jobs.
  -local            Combination of 'localHost' and 'localDisk'
  -localHost        Assemble on this host (default: on grid)
  -localDisk        Use current directory to assemble. (default: aserver)
  -D <name>         Specify project name (displayed on Aserver Console).
                    Use -D test for debugging. (default: frg prefix) 
  -[no]notify    	Send email upon completion/error and create BITS case 
                    for errors. (default: notify)
  -after            Print the generated shell script and exit
  -cont <req_id>    Use the request id given for continuation.
  -invocation <s>   Original user invocation (for logging purposes)
  -service <s>      Provide the service type:
                    CA - for Celera Assemblies
                    TA - for TIGR Assemblies
                    POST - for Post Assembly processes                    
  -test             Run in debug mode (same as -D test -nonotify)  
  -wait             Wait for job completion before exiting.
       
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
use DBI;
use Cwd;               #Use cwd to retrive current working directory.
use File::Basename;    #Use fileparse to get path, name, and extension of file
use Sys::Hostname;     #Use hostname
use TIGR::Foundation;
use TIGR::ConfigFile;
use IO::File;
use File::Spec;        # rel2abs()
use File::Copy;        # copy()

# ================================== Constants  ===============================
# default config file for arun
my $DEFAULT_INSTALL_DIR = "/usr/local/common/CARUN";
my $ARUN_CONFIG_FILE = "ARUN.conf";
my $install_dir = $ENV{'ARUN_INSTALL_DIR'};    

# The arun Config file object
my $arun_cf = undef;

#Config file parameters set by initializeConig
my $jcvi             = undef;
my $debug            = undef;
my $CMD_LOG          = undef;
my $CMD_ERR          = undef;
my $sybase           = undef;
my $server           = undef;
my $type             = undef;
my $user             = undef;
my $pass             = undef;
my $database         = undef;
my $host             = undef;
my $marshalling_base = undef;
my $props_file       = undef;
my $support_email    = undef;
my $config_completed = 0;

# when Cntl-C is used to cancel arun
my $cancel_track = 0;

#Number of steps used by ca_observer
my $cmd_count = 0;

#request id from asdb for this assembly run.
my $request_id = undef;

#aget default modes
my $AGET_MIN_MODE = 2;
my $AGET_MED_MODE = 1;
my $AGET_MAX_MODE = 0;

#default service name
my $SERVICE_CELERA_ASSEMBLY = 'ASSEMBLY';
my $SERVICE_TIGR_ASSEMBLY = 'TIGR_ASSEMBLY';
my $SERVICE_POST_ASSEBMLY = 'POST_ASSEMBLY';
my $DEFAULT_SERVICE_NAME = $SERVICE_CELERA_ASSEMBLY;

my %SERVICE_HASH = (
    CA => $SERVICE_CELERA_ASSEMBLY,
    TA => $SERVICE_TIGR_ASSEMBLY,
    POST => $SERVICE_POST_ASSEBMLY
);

#grid executable name
my $generatedScript = 'command.sh';
my $inputRecipe = 'recipe.sh';
# database handle
my $dbh = undef;

# Module for housekeeping functions
my $tf_object = new TIGR::Foundation;

# User directory from which arun is invoked
my $basedir = $tf_object->getProgramInfo("env_path");

# Clear mask (for future file/directory creation)
umask(0000);

# Required files to copy
my @requiredFiles = ();

# Recipe lines read into memory
my @recipeLines = ();


#Name:   initializeConfig
#Input:  config file name
#Usage:  Retrive config file values
sub initializeConfig($) {

    my $configFile = shift;
    $install_dir = $DEFAULT_INSTALL_DIR if ( !defined $install_dir );    
    $configFile = $install_dir . '/' . $ARUN_CONFIG_FILE if ( !defined $configFile );

    # The arun Config file object
    $arun_cf = new TIGR::ConfigFile($configFile)
      or bail("Could not initialize the config file: '$configFile'");
    $tf_object->logLocal( "Using config file: $configFile", 2 );

    if( $arun_cf->getOption('jcvi') == 1) {
        $jcvi = 1;
        $tf_object->logLocal( "Running at JCVI (GRID available)", 1 );
    } else {
        $tf_object->logLocal( "Running offsite", 1 );
    }

    # logging for the generated executing script
    $CMD_LOG          = $arun_cf->getOption('command_log');
    $CMD_ERR          = $arun_cf->getOption('command_err');
    $marshalling_base = $arun_cf->getOption('marshalling_base');

    if ( $jcvi == 1 ) {
        # Default parameters
        $sybase           = $arun_cf->getOption('sybase','jcvi');
        $server           = $arun_cf->getOption('server','jcvi');
        $type             = $arun_cf->getOption('type','jcvi');
        $user             = $arun_cf->getOption('db_user','jcvi');
        $pass             = $arun_cf->getOption('db_pass','jcvi');
        $database         = $arun_cf->getOption('database','jcvi');
        $host             = $arun_cf->getOption('host','jcvi');
        $props_file       = $marshalling_base .'/'. $arun_cf->getOption("propfile",'jcvi');
        $support_email    = $arun_cf->getOption('support_email','jcvi');
    }
    
    $config_completed = 1;    
}

sub jcviOptions($$$$$$$$) {
    my $alias = shift;
    my $cancel_reqId = shift;
    my $maxDest = shift;
    my $minDest = shift;
    my $medDest = shift;
    my $custCopyStr = shift;
    my $fast = shift;
    my $notify = shift;
        
    bail("Unsupported option ('jcvi' flag must be set).")
        if ( defined $alias or defined $cancel_reqId or defined $maxDest 
                or defined $minDest or defined $medDest or defined $custCopyStr 
                or defined $fast or defined $notify);
}

#Name:   prepare_script
#Input:  See comments per variable
#Usage:  Generate the output script of the input recipe
sub prepare_script($$$$$) {
    my $prefix          = shift;
    my $origdir         = shift; #User calling directory
    my $outdir          = shift; #Working directory under /local/aserver
    my $resdir          = shift; #Aget directory
    my $getMode         = shift; #Indicates level of aget, if any

    #Executables
    my $scriptRunningScript = $install_dir . '/' . $arun_cf->getOption('asdbrunning_script','jcvi');
    my $ca_observer         = $arun_cf->getOption('observer_executable','jcvi');
    my $aget                = $arun_cf->getOption('aget_executable','jcvi');

    my $ca_observerLog = "$outdir/log/ca_observer.log";

    ### ##########################################
    ### Create scripts directory, this is
    ### where the output scripts will go
    ### ##########################################
    my $output_file = "$outdir/scripts/$generatedScript";

    ### ##########################################
    ### Add Header Information, Variables
    ### to output script
    ### ##########################################
    my $whoami         = `whoami`; chomp($whoami);
    my $generationTime =   `date`; chomp($generationTime);
    my $invocation     = $tf_object->getProgramInfo("invocation");

    my $output_fh = new IO::File ">$output_file"    
      or bail("Could not create executable file: '$output_file'. Error Code: $!");
    print $output_fh "#!/bin/sh\n";
    print $output_fh 
"echo \"arun_header> Arun Generated file by $whoami on $generationTime\"\n";
    print $output_fh "echo \"arun_header> Arun Invocation: $0 $invocation\"\n";
    print $output_fh "echo \"arun_header> Arun Version: $MY_VERSION\"\n";
    print $output_fh "echo \"arun_header> User Directory: $origdir\"\n";
    print $output_fh "echo \"arun_header> Work Directory: $outdir\"\n";
    print $output_fh "echo \"arun_header> Running on `hostname`\"\n";
    print $output_fh 
"echo \"arun_header> --------------------------------------------------------\"\n";
    print $output_fh "cp $origdir/arun.log $outdir/log/.\n";
    
    if ( $jcvi == 1) {
        print $output_fh "#Make sure that SYBASE is defined in the environment.\n";
        print $output_fh "echo \$SYBASE | grep -i sybase > /dev/null 2>&1\n";
        print $output_fh "STATUS=\$?\n";
        print $output_fh "if [ \$STATUS != 0 ]; then\nSYBASE=$sybase\nexport SYBASE\n";
        print $output_fh "fi\n";
    }

    #Set variables
    print $output_fh "echo \"arun_header> Environment Variables\"\n";
    print $output_fh "ORIGDIR=$origdir\n";
    print $output_fh "echo ORIGDIR=$origdir\n";
    print $output_fh "WORKDIR=$outdir\n";
    print $output_fh "echo WORKDIR=$outdir\n";
    print $output_fh "PREFIX=$prefix\n";
    print $output_fh "echo PREFIX=$prefix\n";
    print $output_fh "REQUEST_ID=$request_id\n";
    print $output_fh "echo REQUEST_ID=$request_id\n";
    print $output_fh "export ARUN_INSTALL_DIR=$install_dir\n";
    print $output_fh "echo REQUEST_ID=$install_dir\n";

    #Add runScripts to output script to update asdb when the script hits the grid
    unshift( @recipeLines, "$scriptRunningScript $request_id\n" )
        if ( $jcvi == 1);

    #new size of recipe lines to execute
    my $size = @recipeLines;
    
    #Description of the current command, used for logging and error handling
    my $command = 'Initializing';
    for ( my $counter = 0 ; $counter < $size ; $counter++ ) {
        my ($line) = $recipeLines[$counter];
        chomp($line);

        if ( $line =~ /^\s*#>\s*(.+)/ )
        {    #Command Comment (Seen as a 'Step' on the Aserver Console)
            if ( $cmd_count > 0 ) {
                if ( $jcvi == 1 ) {
                    print $output_fh "$ca_observer --appendlog=1 --logfile=$ca_observerLog --event=finish --name=\"$command\" --retval=0 --props=$props_file --host=`hostname` --message=\"Command with name: '$command' finished\"\n"
                } else {
                    print $output_fh "echo \"Step $cmd_count: '$command' finished\"";
                }                
            }
            $cmd_count++;
            $command = $1;
            print $output_fh  "echo \"arun_step> $command \"\n";
            print $output_fh "$ca_observer --appendlog=1 --logfile=$ca_observerLog --event=start --name=\"$command\" --retval=0 --props=$props_file --host=`hostname` --message=\"Command with name: '$command' started\"\n"
                if ( $jcvi == 1 );
        }
        elsif ( $line =~ /^\s*exit/ )
        {    #If exit in recipe it will stop reading the recipe file
            $counter = $size;
        }
        elsif ( $line =~ /^#!\/bin\/sh\s*$/ ) {    #ignores the shebang line
            next;
        }
        elsif ( $line =~ /^\s*#\s*(.+)$/ ) {    #Prints out Comment lines
            $command = $1;
            print $output_fh  "echo \"arun_cmt> $command \"\n";
        }
        elsif ( $line =~ /^\s*$/ ) {            #If blank line, skip it
            next;
        }
        else {  #Else it is an executable line
            #Check if command spans multiple lines using '\' or '&&'
            #If so, combines them into one line
            while ( $line =~ /^(.*)(?:(&&)|\\)\s*$/ ) {
                $line = $1;
                if ( defined $2 ) {
                    $line .= ' ' . $2 . ' ';
                }
                $counter++;
                if ( $counter == $size ) {
                    last;
                }
                my $tmpLine = $recipeLines[$counter];
                chomp $tmpLine;
                $line .= ' ' . $tmpLine;
            }

            my $lineEcho = $line;
            $lineEcho =~ s/\$/\\\$/g;
            print $output_fh "echo \"arun_cmd> '$lineEcho'\"\n";
            print $output_fh "TIME_INIT=`/bin/date +\%s`\n";
            print $output_fh "$line\n";
            print $output_fh "STATUS=\$?\n";
            print $output_fh "echo \"arun_stat> STATUS=\$STATUS\"\n";
            print $output_fh "TIME_FINAL=`/bin/date +\%s`\n";
            print $output_fh 
              "EXECTIME=`/usr/bin/expr \$TIME_FINAL - \$TIME_INIT`\n";
            print $output_fh 
              "echo \"arun_time> Done. Took \$EXECTIME seconds.\"\n";
            print $output_fh 
"if [ \$STATUS != 0 ]; then\necho \"arun_stat> ERROR, exit status = \$STATUS\"\n";
            print $output_fh "$ca_observer --appendlog=1 --logfile=$ca_observerLog  --event=failure --name=\"$command\" --retval=0 --props=$props_file --host=`hostname` --message=\"Command with name: '$command' failed\"\n"
                if ( $jcvi == 1 );
            print $output_fh "STATUS=$?\n";
            print $output_fh "chmod 664 $outdir/$CMD_LOG\n";
            print $output_fh "cp $outdir/log/command.* $resdir\n" if ( $outdir != $resdir );
            print $output_fh "exit\n";
            print $output_fh "fi\n";
        }
    }

    print $output_fh "$ca_observer --appendlog=1 --logfile=$ca_observerLog --event=finish --name=\"$command\" --retval=0 --props=$props_file --host=`hostname` --message=\"Command with name: '$command' finished\"\n"
        if ( $jcvi == 1 );

    #If getMode set, call to aget will be added to the recipe (but not counted as one of the steps)
    if ( defined($getMode) and $jcvi == 1 ) {
        my $agetCmd;
        if ( $getMode ne $AGET_MAX_MODE and 
             $getMode ne $AGET_MED_MODE and 
             $getMode ne $AGET_MIN_MODE ) {     #CUSTOMIZE, $getMode is the list file name
           $agetCmd =  "$aget -L $getMode";
        } else {
           $agetCmd =  "$aget -mode $getMode";
        }        
        $agetCmd .= "  $request_id $resdir >> $outdir/$CMD_LOG";
        
        print $output_fh  "cd $outdir\n";
        print $output_fh "$ca_observer --appendlog=1 --logfile=$ca_observerLog --event=start --name=aget --retval=0 --props=$props_file --host=`hostname` --message=\"Command with name: 'aget' started\"\n";
        print $output_fh "$agetCmd\n";
        print $output_fh "STATUS=\$?\n";
        print $output_fh "echo \"arun_stat> STATUS=\$STATUS\"\n";
        print $output_fh "if [ \$STATUS != 0 ]; then\necho \"arun_stat> ERROR, exit status = \$STATUS\"\n";
        print $output_fh "$ca_observer --appendlog=1 --logfile=$ca_observerLog  --event=failure --name=aget --retval=0 --props=$props_file --host=`hostname` --message=\"Command with name: aget failed\"\n";
        print $output_fh "STATUS=$?\n";
        print $output_fh "chmod 664 $outdir/$CMD_LOG\n";
        print $output_fh "cp $outdir/log/command.* $resdir\n";
        print $output_fh "exit\n";
        print $output_fh "fi\n";
    }
    print $output_fh  "chmod 664 $outdir/$CMD_LOG\n";
    print $output_fh  "exit\n";
    close $output_fh 
      or bail("Error closing '$output_file'. Error Code: $!");
    chmod 0755, $output_file
        or die ("Unable to change mode: $!");
    return $output_file;
}

#Name:   init_prop_file
#Usage:  Initialize the properties file
sub init_prop_file($$$) {
    my $resdir    = shift;
    my $notify    = shift;
    my $prefix    = shift;

    my $props_fh = new IO::File ">$props_file"
      or bail("Could not open props file: '$props_file'. Error Code: $!");
    print $props_fh qq~
    request_id=$request_id
    server=$server
    user=$user
    password=$pass
    database=$database
    host=$host
    marshalling_base=$marshalling_base
    clean=9
    prefix=$prefix
    total_commands=$cmd_count
    command_count=0
    ~;
    if ( defined $resdir ) {
        print $props_fh "resdir=$resdir\n";
    }
    if ( $notify == 1 ) {
        print $props_fh "email=$notify\n";
        print $props_fh "support_email=$support_email\n";
    }
    close $props_fh 
      or bail("Error closing '$props_file'. Error Code: $!");
}

#Name:   init_db_connection
#Input:  none
#Output: none
#Usage:  This function initializes a connection to the asdb database.
sub init_db_connection() {
    $tf_object->logLocal( "Establishing a connection to the database", 3 );

    # Try to connect to the database server
    $dbh = DBI->connect(
        "dbi:$type:server=$server;packetSize=8092",
        $user, $pass,
        {
            PrintError => 0,
            RaiseError => 0,
            AutoCommit => 1
        }
    );
    $dbh->do("use $database")
      or $tf_object->bail("Failed to open database \'$database\'");
    $dbh->{InactiveDestroy} = 1;
}

#Name:   close_db_connection
#Input:  none
#Output: none
#Usage:  This function initializes a connection to the asdb database.
sub close_db_connection() {
    $tf_object->logLocal( "Closing connection to the database", 3 );
    $dbh->disconnect();
    $dbh = undef;
}

#Name:   asdbInit
#Input:  project name and alias for request id
#Output: none
#Usage:  ASDB initializing for the new request
sub asdbInit ($$$) {

    my $project = shift;
    my $alias   = shift;
    my $service = shift;

    #First,
    #select service_type_id from Service_Type

    init_db_connection() unless defined $dbh;

    my $query =
"select service_type_id from Service_Type where service_type = '$service'";

    $tf_object->logLocal( "Executing the query $query", 3 );

    my $qh = $dbh->prepare($query)
      or bail("Cannot prepare $query: " . $dbh->errstr );

    defined $qh->execute()
      or bail("Database query \'$query\' failed: " . $dbh->errstr );

    my @row             = $qh->fetchrow();
    my $service_type_id = $row[0];
    $tf_object->logLocal( "Service_type_id=$service_type_id", 3 );
    $qh->finish();
    $qh = undef;

    #Second
    #insert into Request

    $query =
"insert into Request (status,job_entered,submitter,submitter_host,service_type_id,job_started,job_terminated,submission_type,alias) "
      . "values ('I',getdate(),'"
      . getpwuid($<) . "','"
      . hostname()
      . "',$service_type_id,"
      . 'convert(smalldatetime,getdate()),' . 'NULL,' . '1';

    if ( defined $alias ) {
        $query .= ",'$alias')";
    }
    else {
        $query .= ",NULL)";
    }

    $tf_object->logLocal( "Running query $query ...", 2 );
    my $qh = $dbh->prepare($query)
      or bail("Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or bail("Database query \'$query\' failed: " . $dbh->errstr );

    $qh->finish();
    $qh = undef;

    #Third
    #get generated request id
    $query = "select \@\@identity";

    $tf_object->logLocal( "Executing the query $query", 3 );

    my $qh = $dbh->prepare($query)
      or bail("Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or bail("Database query \'$query\' failed: " . $dbh->errstr );

    my @row        = $qh->fetchrow();
    my $request_id = $row[0];
    $tf_object->logLocal( "request_id=$request_id", 3 );
    $qh->finish();
    $qh = undef;
    $tf_object->logLocal( "The aserver request_id for this job is $request_id",0 );

    #Fourth
    #update project name in asdb
    if ( $project eq 'test' ) {
        $project = "test\n";
    }
    $query =
        "insert into Project (request_id, project_name) "
      . "values ($request_id, \"$project\")";

    $tf_object->logLocal( "Running query $query ...", 2 );
    my $qh = $dbh->prepare($query)
      or bail( "Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or bail( "Database query \'$query\' failed: " . $dbh->errstr );

    $qh->finish();
    $qh = undef;

    #Fifth
    #Make directory request directory
    my $request_dir = $marshalling_base . '/' . $request_id;
    mkdir "$request_dir", 0777 unless -d "$request_dir";
    $tf_object->logLocal( "request_dir=$request_dir", 3 );

    #Sixth
    #Update progress in Process table
    $query =
        "insert into Process (request_id,pid,progress,cluster_id) "
      . "values ($request_id,1,0,0)";

    $tf_object->logLocal( "Running query $query ...", 2 );
    my $qh = $dbh->prepare($query)
      or bail( "Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or bail( "Database query \'$query\' failed: " . $dbh->errstr );

    $qh->finish();
    $qh = undef;

    return $request_id;
}

#Name:   printRecipeAfter
#Input:  path to generated script
#Output: none
#Usage:  Print generated script to STDOUT
sub printRecipeAfter($) {
    my $script = shift;

    my $script_fh = new IO::File "<$script"
      or bail( "Could not open script file: '$script'. Error Code: $!" );

    my @script_array = <$script_fh>; 
    foreach (@script_array) {
        print STDOUT $_;
    }    
    close $script_fh 
      or bail("Error closing '$script'. Error Code: $!");
    
    setStatusFinished();
    exit 0;
}


#Name:   setStatusFinished
#Input:  none
#Output: none
#Usage:  Indicate to ASDB that request is complete
sub setStatusFinished () {
    if ( $config_completed == 1 and defined $request_id and $jcvi == 1 ) {
        init_db_connection() unless defined $dbh;
    
        #Set Status
        my $start_query =
    "update Request set status = \'F\', job_terminated = convert(smalldatetime,getdate()) where request_id = $request_id";
        $tf_object->logLocal( "Executing the status query $start_query", 3 );
        my $qh = $dbh->prepare($start_query)
          or bail( "Cannot prepare $start_query: " . $dbh->errstr );
        defined $qh->execute()
          or bail( "Database query \'$start_query\' failed: " . $dbh->errstr );
    
        #Set Progress
        my $progress_query =
          "update Process set progress = 100 where request_id = $request_id";
        $tf_object->logLocal( "Executing the progress query $progress_query", 3 );
        my $qh = $dbh->prepare($progress_query)
          or
          bail( "Cannot prepare $progress_query: " . $dbh->errstr );
        defined $qh->execute()
          or bail( "Database query \'$progress_query\' failed: " . $dbh->errstr );
    }
}

#Name:   printInvocationScript
#Input:  none
#Output: none
#Usage:  Create invocation script
sub createInvocationScript($) {
    my $invocation = shift; 
    my $invocation_script = $arun_cf->getOption('invocation_script');
    #writing the invocation
    my $invo_fh = new IO::File ">$invocation_script"
      or bail("Cannot open '$invocation_script'. Error code: $!");
     
    if ( !defined $invocation or $invocation =~ /^\s*$/ ) {
        $invocation = $0,$tf_object->getProgramInfo("invocation");
    }
    
    my $cust_name  = getpwuid($<);
    my $hostname   = hostname();
    print $invo_fh "invocation: $invocation\n";
    print $invo_fh "username: $cust_name\n";
    print $invo_fh "hostname: $hostname\n";
    print $invo_fh "userdir: $basedir\n";
    close $invo_fh 
      or bail("Error closing '$invocation_script'. Error Code: $!");

}

#Name:   readRecipe
#Input:  recipe_file name
#Output: none
#Usage:  Reads recipe file into memory 
sub readRecipe() {    
    @recipeLines = <STDIN>;
}

#Name:   writeRecipeToFile
#Input:  working_dir
#Output: none
#Usage:  Output the input recipe to file (for logging purposes) 
sub writeRecipeToFile($) {
    my $working_dir = shift;
    
    my $recipeFilePath = "$working_dir/scripts/$inputRecipe";
    my $recipe_fh = new IO::File (">$recipeFilePath")
      or bail("Cannot write to: '$recipeFilePath'. Error code: $!");
    foreach (@recipeLines) {
        print $recipe_fh $_;
    }
    close $recipe_fh;
}
        


#Name:   checkRequiredFiles
#Input:  user_dir, prefix
#Output: none
#Usage:  Checks the required files from the recipe 
sub checkRequiredFiles($$) {
    
    my $user_dir = shift;
    my $prefix = shift;
    
    foreach my $recipeLine (@recipeLines) {
        if ( $recipeLine =~ /^##REQUIRED\s*(.+)/ ) {
            my $requiredStr = $1;
            $requiredStr =~ s/PREFIX/$prefix/g;
            my @requiredArray = split(/\s+/, $requiredStr);
            foreach my $requiredFile (@requiredArray) {
                if ( $requiredFile !~ /^\s*$/ ) {                    
                    if ( !-r "$user_dir/$requiredFile") {
                        my $tmpFile = $requiredFile;
                        $tmpFile =~ s/$prefix/tmp/g;
                        my $bailMsg = "Missing required file '$user_dir/$requiredFile' for specified recipe."; 
                        if ( -r "$user_dir/$tmpFile") {
                            $bailMsg .= "\nBut found $tmpFile.  Did you mean to run with 'tmp' instead of '$prefix'?";
                        }
                        bail($bailMsg); 
                    }
                    push @requiredFiles, $requiredFile;
                }
            }
        }
    }
}

#Name:   copyRequiredFiles
#Input:  user_dir - input directory, working_dir - output directory
#Output: none
#Usage:  Copies the required Files (as determined in checkRequiredFiles) to the output dir 
sub copyRequiredFiles ($$) {
    my $user_dir = shift;
    my $working_dir = shift;    
    foreach my $file (@requiredFiles) {
        copy("$user_dir/$file","$working_dir/$file");
    }    
}

#Name:   prefixProjectSetup
#Input:  project name (from command line)
#Output: user_dir - the directory of the frg file
#        project_out - the project name to use
#        prefix - the prefix of the project to use 
#Usage:  Determines the prefix, project names
sub prefixProjectSetup($) {
    my $project_in = shift;
    
    my $prefix = undef;
    my $project_out = undef;
    my $user_dir = undef;
    
    # Setup assembly prefix
    my $argvSize = scalar(@ARGV);
    
    if ( $argvSize != 1 ) {
        bail("Please specify an assembly prefix.");
    }
    elsif ( $argvSize == 1 ) {
        my $fileentry = $ARGV[0];

        if ( $fileentry !~ /\.frg$/ ) {
            $fileentry .= ".frg";
        }
        $tf_object->logLocal( "Single frg: $fileentry", 2 );

        my ( $name, $path, $suffix ) = fileparse( $fileentry, ".frg" );
        $user_dir = File::Spec->rel2abs($path);
        $prefix = $name;
    }

    # Set up for project name to be sent to asdb
    # If not defined via command-line, get from .project file,
    # else get from ASSEMBLY Prefix
    if ( !defined $project_in ) {
        my $project_file = "$prefix.project"; 
        if ( -e "$project_file" ) {
            my $project_fh = new IO::File "<$project_file"
                or bail("Cannot open: '$project_file'. Error code: $!");
            my @prjFileLines = <$project_fh>;
            chomp( $project_out = $prjFileLines[0] );
            close $project_fh 
                or bail("Error closing '$project_file'. Error Code: $!");
        }
        else {
            $project_out = $prefix;
        }
    }
    else {
        $project_out = $project_in;    
    }

    return $user_dir, $project_out, $prefix;
}

#Name:   setupWorkDirectory
#Input:  request_id - request id of current job 
#Output: working_dir - path to the working directory
#Usage:  Setup the working directory, copies over continue directory if provided 
sub setupWorkDirectory($$$$$) {

    my $request_id = shift;
    my $user_dir = shift;
    my $continue_ReqId = shift;
    my $localDisk = shift;
    my $curr_dir = shift; #Typically is the same as user_dir, but in case user is running
                          #from a different directory than that where the original files are located
                          #(user_dir)
        
    my $working_dir = $localDisk ? $curr_dir :"$marshalling_base/$request_id";
        
    my $progress = "Assembly output being written to \'$working_dir\'";
    print STDOUT "$progress...\n";
    $tf_object->logLocal( "$progress", 1 );

    #Update props_file -> later used by init_prop_file, prepare_script
    $props_file =~ s/REQUEST_ID/$request_id/g;
    $tf_object->logLocal( "Using props file: $props_file", 2 );


    $CMD_LOG =~ s/REQUEST_ID/$request_id/g;
    $CMD_ERR =~ s/REQUEST_ID/$request_id/g;
    
    if (!defined $continue_ReqId) {
        # Copy the required (unless continuation)
        copyRequiredFiles($user_dir,$working_dir) if ( !$localDisk );    
        
        # creating the log directory for storing the log files
        mkdir( "$working_dir/log", 0777 ) unless -d "$working_dir/log";
        # create scripts directory
        mkdir "$working_dir/scripts", 0777 unless -d "$working_dir/scripts";
        
    } else {
        $tf_object->runCommand("cp -r $marshalling_base/$continue_ReqId/* $marshalling_base/$request_id/.");
        unlink $props_file,
               "$marshalling_base/$request_id/scripts/$generatedScript";
    }
    
    return $working_dir;
}

#Name:   continueValidation
#Input:  continue_ReqId - request if job to continue 
#Output: none
#Usage:  Validate the continuation directory 
sub continueValidation($) {
    my $continue_ReqId = shift;
    
    bail ("Continuation not possible '$marshalling_base/$continue_ReqId' does not exist.")
        if ( !-e "$marshalling_base/$continue_ReqId");
}

#Name:   cdToWorkingDirectory
#Input:  working_dir
#Output: none
#Usage:  Change directory to working directory 
sub cdToWorkingDirectory($) {
    my $working_dir = shift;
    
    # Change directory to working directory
    chdir $working_dir or bail("Failed to cd to $working_dir");
    $tf_object->logLocal( "changing to \'$working_dir\'", 1 );
}


#Name:   validateServiceType
#Input:  $service - Service Name
#Output: none
#Usage:  Validates the service type option and sets it to default if not set 
sub validateServiceType($) {
    my $service = shift;
    
    if ( !defined $service) {
        $service = $DEFAULT_SERVICE_NAME;
        return;
    } elsif ( !exists $SERVICE_HASH{$service}) {
        bail("Invalid service type: '$service'");
    } else {
        $service = $SERVICE_HASH{$service};
    }
    
    return $service;
}


#Name:   checkAgetOptions
#Input:  Copy options
#Output: none
#Usage:  Checks for violations of aget options (which are mutually exclusive)). 
sub checkAgetOptions($$$$$$) {
    my $noCopy = shift;
    my $minDest = shift;
    my $medDest = shift;
    my $maxDest = shift;
    my $custCopyStr = shift;
    my $localDisk = shift;

    if ( defined $localDisk and 
        (defined $minDest or defined $medDest or defined $maxDest or defined $custCopyStr) ) {
        bail("'-localDisk' option is not compatible with '---Copy' options");
    }
    my $numDefined = 0;
    $numDefined++ if( defined $noCopy ); 
    $numDefined++ if( defined $minDest ); 
    $numDefined++ if( defined $maxDest );
    $numDefined++ if( defined $medDest );
    $numDefined++ if( defined $custCopyStr );
    
    bail("The options: 'noCopy' 'minCopy' 'medCopy' 'maxCopy' 'custCopy' are mutually exclusive, please choose one.")
        if ( $numDefined > 1 );
}    

#Name:   agetSetup
#Input:  minDest - command-line option
#        medDest - command-line option
#        maxDest - command-line option
#Output: getMode - aget mode
#        copyDir - directory to copy to
#Usage:  Setup variables for calling aget 
sub agetSetup($$$$$) {
    my $user_dir    = shift;
    my $minDest     = shift;
    my $medDest     = shift;
    my $maxDest     = shift;
    my $custCopyStr = shift;

    #Default output directory is the user invocation directory
    my $copy_dir = File::Spec->rel2abs($user_dir);
    my $getMode = $AGET_MED_MODE;
       
    if ( defined $custCopyStr ) {
        my ($custCopyFile,$custCopyDir) = split (',',,$custCopyStr);
        bail ("Unable to read file: $custCopyFile") if ( !-r $custCopyFile );
        $copy_dir = File::Spec->rel2abs($custCopyDir)  if ( defined $custCopyDir and $custCopyDir ne '' );
        $getMode = $custCopyFile; 
    } 
    elsif ( defined $minDest ) {
        $copy_dir = File::Spec->rel2abs($minDest)  if ( $minDest ne '' );
        $getMode = $AGET_MIN_MODE;
    }
    elsif ( defined $maxDest ) {
        $copy_dir = File::Spec->rel2abs($maxDest) if ( $maxDest ne '' );
        $getMode = $AGET_MAX_MODE;
    }
    elsif ( defined $medDest ) {
        $copy_dir = File::Spec->rel2abs($medDest) if ( $medDest ne '' );
    }  

    # Create local directory for aget, if specified
    if ( defined($getMode) && ( !-e $copy_dir ) ) {
        mkdir( $copy_dir, 0777 ) unless -d "$copy_dir";
    }
    
    return $getMode, $copy_dir;
}

#Name:   setDebugLevel
#Input:  none
#Output: none
#Usage:  Use command-line option and config value to set debug level 
sub setDebugLevel() {
    # Update debug level for logging if it is not specified by the user
    $tf_object->setDebugLevel($arun_cf->getOption('debug_default'))
        if ( !defined $tf_object->getDebugLevel() );             
    $debug = $tf_object->getDebugLevel();
}

#Name:   executeScript
#Input:  offgrid - specifies to run locally 
#        shell_script - script to run
#        working_dir - directory where script will run
#        getMode - aget mode
#        copy_dir - directory for aget to copy to
#        wait - wait for job completion before exiting
#Output: none
#Usage:  Executing the generated script 
sub executeScript($$$$$$$) {

    my $offgrid = shift;
    my $shell_script = shift;
    my $working_dir = shift;
    my $getMode = shift;
    my $copy_dir = shift;
    my $wait = shift;
    my $fast = shift;
    
    my $ca_cmd;
    my $runMsg = "Running the command \'$ca_cmd\'";
    my $msg    = "Running the CA request";
    my $runCommandLogDir;

    if ( $offgrid ) {
        $ca_cmd = $shell_script;
        $runMsg .= ' locally.';
        $msg    .= " locally.\n";
        $runCommandLogDir = "$working_dir/$CMD_LOG";
    }
    else {
        $ca_cmd =
          $arun_cf->getOption('grid_cmd','jcvi')
          . " -o $working_dir/$CMD_LOG -e $working_dir/$CMD_ERR";
        $ca_cmd .= ' ' . $arun_cf->getOption('grid_sync','jcvi') . ' ' if ($wait);
        $ca_cmd .= ' ' . $arun_cf->getOption('grid_fast','jcvi') . ' ' if ($fast);        
        $ca_cmd .= " $shell_script";
        $runMsg .= ' on the grid.';
        $msg    .= " on the grid.\n";
        $runCommandLogDir = "$working_dir/" . $arun_cf->getOption('submission_log') ;
    }

    if ( defined($getMode) ) {
        my $copy_dir_name = basename($copy_dir);
        $msg .= "The local output directory is $copy_dir_name.\n";
    }
    $msg .= "The request id is $request_id.\n";
    print STDOUT $msg;

    my $loggingRedirects = " >> $runCommandLogDir 2> $working_dir/$CMD_ERR";
    $ca_cmd .= $loggingRedirects;
    $tf_object->runCommand("echo \"\$runMsg\" $loggingRedirects");
    $tf_object->logLocal( "$runMsg", 1 );
    $tf_object->logLocal( "$ca_cmd", 1 );
    my $bad = $tf_object->runCommand("$ca_cmd");
    if ($bad) {
        bail("The command $ca_cmd was not successful");
    }
    elsif ( $offgrid ) {
        print STDOUT
"Your job has completed.\n";
    }
    else {    #Successful submission to grid
        print STDOUT
          "You can cancel the job by running: arun -cancel $request_id\n";
    }
    print STDOUT
"Please check the console for your job's status at \nhttp://aserver.tigr.org:8080/AserverConsole/\n";
    
}

# =============================== MAIN ======================================
#
MAIN:
{
    my $alias           = undef;    # Project alias in ASDB
    my $after           = undef;    # Option to generate and print but not execute script
    my $local           = undef;    # Combine offgrid and localDisk
    my $offgrid         = undef;    # Run locally
    my $localDisk       = undef;    # Use current directory to execute
    my $version         = undef;    # Version option
    my $prefix          = undef;    # Prefix of required files
    my $project         = undef;    # Specify project name for ASDB 
    my $maxDest         = undef;    # maxCopy option (optionally dir name)    
    my $minDest         = undef;    # minCopy option (optionally dir name)
    my $medDest         = undef;    # medCopy option (optionally dir name)
    my $noCopy          = undef;    # noCopy option
    my $custCopyStr     = undef;
    my $custCopyFile    = undef;
    my $custCopyDir     = undef;
    my $cancel_reqId    = undef;    # cancel job option
    my $getMode         = undef;    # agetMode    
    my $configFile      = undef;    # optional user specified config file
    my $continue_ReqId  = undef;
    my $invocation      = undef;    # User invocation of calling program
    my $service         = undef;             
    my $notify          = 1;        # Enabled by default
    my $test            = undef;
    my $wait            = undef;    # Wait for job completion before exiting
    my $fast            = undef;    
    
    # Directory variables
    my $working_dir = undef;    #location the recipe file will execute out of
    my $copy_dir = undef;       #location where aget will copy files
    my $curr_dir = cwd;
    my $user_dir = undef;         #location from which arun was invoked

    # ========================== Program Setup ==============================
    # Prepare logs
    $tf_object->addDependInfo(@MY_DEPENDS);
    $tf_object->setHelpInfo($HELPTEXT);
    $tf_object->setVersionInfo($MY_VERSION);

    $tf_object->TIGR_GetOptions(
        'alias=s',   \$alias,           'cancel=i',  \$cancel_reqId,
        'minCopy:s', \$minDest,         'maxCopy:s', \$maxDest,
        'medCopy:s', \$medDest,         'noCopy',    \$noCopy,
        'custCopy=s',\$custCopyStr,     'local',     \$local,
        'offgrid|localHost', \$offgrid, 'localDisk', \$localDisk,  
        'after',      \$after,          'fast',      \$fast,
        'notify!',   \$notify,          'version|V', \$version,
        'D=s',       \$project,         'config=s',  \$configFile,
        'cont=i',    \$continue_ReqId,  'invocation=s', \$invocation,
        'test',      \$test,            'wait', \$wait,
        'service=s', \$service
    ) or bail("Options could not be read. See arun -h.");

    # Initialize parameters based on config file
    initializeConfig($configFile);

#    jcviOptions($alias,$cancel_reqId,$maxDest,$minDest,$medDest,$custCopyStr,$fast,$notify)
#        if ( !$jcvi );

    if ( $local or !$jcvi ) {
        $offgrid = 1;
        $localDisk =1 ;
    }

    # Check aget options - ensure mutual exclusivity
    checkAgetOptions($noCopy,$minDest,$medDest,$maxDest,$custCopyStr,$localDisk)
        if ( $jcvi == 1 );
    
    #validate service type
    $service = validateServiceType($service)
        if ( $jcvi == 1 );
    # Set debug level
    setDebugLevel();
    
    # User request to cancel job for a request_id
    cancelJob($cancel_reqId) if (defined $cancel_reqId && $jcvi == 1);

    # If continuation requested, validate the given request id
    continueValidation($continue_ReqId) if(defined $continue_ReqId && $jcvi == 1);

    if ( defined $test) {
        $notify = 0;
        $project = 'test';
    }  
    # Set up prefix and $project 
    ($user_dir, $project, $prefix) = prefixProjectSetup($project);
    
    # When testing, set the alias to the project name (unless an alias is already provided),
    # this helps to identify the test projects on the aserver console, otherwise, they all
    # show up as 'test'
    $alias = $prefix if (defined $test and !defined $alias);

    # Read recipe into memory
    readRecipe();
        
    # Check that required recipe files exist
    checkRequiredFiles($user_dir,$prefix);    
    
    # Aget setup, unless noCopy specified
    ($getMode, $copy_dir) = agetSetup($user_dir,$minDest, $medDest, $maxDest,$custCopyStr)
        if (!defined $noCopy and !defined $localDisk and $jcvi == 1 );
    
    # Request id setup
    if ($jcvi == 1) {
        $request_id = asdbInit( $project, $alias, $service)
    } else {
        $request_id = time();
    }
    
    # Setup working directory
    $working_dir = setupWorkDirectory($request_id, $user_dir, $continue_ReqId, $localDisk, $curr_dir);
    
    # Change cwd to working directory
    cdToWorkingDirectory($working_dir);
    
    # Create invocation script, used be Aserver Console
    createInvocationScript($invocation);
    
    # Output the input recipe to file (for logging purposes)
    writeRecipeToFile($working_dir);
        
    # Generate the outputted script from the recipe
    my $shell_script = prepare_script(
        $prefix,          $user_dir,     $working_dir,   
        $copy_dir,        $getMode,      
    );

    # If after, generate but do not execute the script
    printRecipeAfter($shell_script) if ($after);
                    
    # Prepare prop file for ca_observer        
    init_prop_file( $copy_dir, $notify, $prefix ) if ($jcvi == 1);

    # Execute script
    executeScript($offgrid,$shell_script,$working_dir,$getMode,$copy_dir,$wait,$fast);

    exit 0;
}

END {
    close_db_connection() if (defined $dbh);
    cancel_request() if ( $cancel_track == 1 );
}

#Cancels a request as long as request is still on the grid.
#Updates asdb to set the status to cancelled, and performs a qdel by finding the
#grid job id in the logs.
sub cancel_request() {
    $tf_object->logLocal( "Call to cancel_request: $config_completed", 3 );

    #If called before config setup is completed, it will not succeed
    if ( $config_completed == 1 and defined $request_id and $jcvi == 1) {
        init_db_connection() unless defined $dbh;

        #Updating the Request table with status "cancelled" for cancelled job
        my $cancel_query = "update Request set job_terminated = convert(smalldatetime,getdate()), "
          . "status = \'C\' where request_id = $request_id";

        $tf_object->logLocal( "Executing the status query $cancel_query", 3 );

        my $qh = $dbh->prepare($cancel_query)
          or $tf_object->bail( "Cannot prepare $cancel_query: " . $dbh->errstr );
        defined $qh->execute()
          or $tf_object->bail( "Database query \'$cancel_query\' failed: " . $dbh->errstr );

        $qh->finish();
        $qh = undef;
    }
}

#Name:   bail
#Input:  Error message
#Output: none
#Usage:  logs error message, calls cancel_request to cancel in ASDB
sub bail ($) {
    my $bailMsg = shift;
    cancel_request();
    $tf_object->bail($bailMsg);
}

#Name:   SIGHANDLER
#Input:  none
#Output: none
#Usage:  this functions acts as the signal handler for the external signals
#        that caexec might receive. It indicates that an external signal has
#        been received by changing the status field in the Request table to
#        "C" (cancel).
sub SIGHANDLER {
    $tf_object->logLocal( "arun received a cancel request", 1 );
    print STDOUT "arun received a cancel request\n";
    $cancel_track = 1;
    exit(2);
}

#Name:   cancelJob
#Input:  cancelId
#Output: none
#Usage:  Cancels a running assembly job on the grid
sub cancelJob ($) {
    my $cancelId = shift;

    $tf_object->logLocal( "Call to cancelJob for jobid: $cancelId", 3 );
    
    my $gridSubmissionLog = "$marshalling_base/$cancelId/" . $arun_cf->getOption('submission_log');    
    
    bail("Couldn't cancel $cancelId, $marshalling_base/$cancelId does not exist.\n")
        if ( !-e "$marshalling_base/$cancelId");
    bail("Couldn't cancel $cancelId, $gridSubmissionLog does not exist.\n")
        if ( !-e "$gridSubmissionLog");
      
    #Retrieved grid job number from log file for this request_id
    my $grepLine =
      `grep \"Your job\" $gridSubmissionLog`;
    chomp($grepLine);
    bail("Couldn't cancel $cancelId, couldn't find qsub ID.\n")
        if ( $grepLine !~ /^Your job (\d+)/ );
        
    my $qdel_result = `qdel $1`;
    $tf_object->logLocal( "Command to cancel: $qdel_result", 3 );    
    chomp($qdel_result);
    bail("The command qdel $1 was not successful")
        if ( $qdel_result =~ /^denied/ );
    print STDOUT "Job $cancelId has been cancelled.\n";
    $request_id = $cancelId;
    $tf_object->logLocal( "Call to cancel_request", 3 );    
    
    cancel_request();
    exit;
}
