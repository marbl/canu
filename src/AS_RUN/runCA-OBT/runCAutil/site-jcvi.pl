use DBI;
use File::Copy;
use FindBin qw($Bin);

my $ca_observer = "$Bin/ca_observer.plx";
my $ca_log      = "log/ca_observer.log";

# These are used only here -- we prefix with local_ just to not
# pollute the namespace too much.

my $local_alias       = undef;
my $local_props_file  = 'props.file';
my $local_maxCopy     = 0;
my $local_medCopy     = 1;
my $local_minCopy     = 2;
my $local_noCopy      = 3;
my $local_copy        = $medCopy;
my $local_notify      = 1;
my $local_test        = 0;

$HELPTEXT =
  qq~Request a whole-genome shotgun assembly using Celera Assembler.

 usage: carun [options] <frg>
 
 inputs: 
  <frg>       Sequence file (.frg) generated e.g. by pullfrag.

 general options:
  -d <dir>          Use <dir> as the working directory. 
                    (default: /usr/local/aserver/var/assembly/<...> )
  -p <prefix>       Use <prefix> as the output prefix. (default: 'asm')
  -alias <a>        Identify the request as <a> on Assembly Server Console
  -maxCopy <dir>    Copy max output to <dir> (default: dir=maxCopy)
  -medCopy <dir>    Copy med output to <dir> (default: dir=medCopy) 
  -minCopy <dir>    Copy min output to <dir> (default: dir=minCopy)
  -noCopy           Do not make local copy. (default: medCopy)
  -[no]notify       Send email upon completion/error and create BITS case 
                    for errors. (default: notify)
  -s <specfile>     Read options from the specifications file <specfile>.
  		    <specfile> can also be one of the following key words:
		    [no]OBT - run with[out] OBT
		    noVec - run with OBT but without Vector
  -test             Run in debug mode (same as -D test -nonotify)  
  -version          Outputs the Celera Assembler Version
  -fields           Outputs the specfile fields
  -help	            This text
 
 CA options: 
  -e <percent>      Unitigger assumed percent error [0.00,0.06] (default: 0.015)
  -g <len>          Assign genome length for unitigger  
  -j <lbound>       Set scaffolder A-stat low bound (default: 1)
  -k <hbound>       Set scaffolder A-stat high bound (default: 5)
  -[no]ubs          Enable unitigger bubble smoothing (default: enabled)  
      
 Genomic Assembly Pipeline: 
 https://intranet.jcvi.org/cms/SE/GAP
 Tracking assembly requests:  
 http://assemblyconsole.tigr.org/
~;



sub localDefaults () {
    
    $wrk = "/usr/local/aserver_new/var/assembly/$request_id" if (!defined($wrk));

    system("mkdir -p $wrk/log") if (! -d "$wrk/log");
    chmod 0755, "$wrk/log";
    caFailure("ERROR: Unable to create log directory.\n") if (! -d "$wrk/log");

    setGlobal("specFile", 'OBT');

    setGlobal("scratch", "$wrk/scratch");

    setGlobal("sge", ' -P 08010 -b n -l msc');
    setGlobal("sgeOverlap", ' -pe threaded 2');
    setGlobal("fakeUIDs", 1);
    setGlobal("useGrid", 1);
    setGlobal("scriptOnGrid", 1);
}



sub localOption($@) {
    my $arg   = shift @_;
    my @ARGV  = @_;
    my $found = 1;

    if ($arg =~ m/^-alias/ ) {
        $local_alias = shift @ARGV;
    } elsif ($arg =~ m/^-maxCopy/ ) {
        $local_copy = $local_maxCopy;
    } elsif ($arg =~ m/^-minCopy/ ) {
        $local_copy = $local_minCopy;
    } elsif ($arg =~ m/^-medCopy/ ) {
        $local_copy = $local_medCopy;
    } elsif ($arg =~ m/^-noCopy/ ) {
        $local_copy = $local_noCopy;
    } elsif ($arg =~ m/^-(no)?notify/ ) {
        if ( $1 eq 'no' ) {
            $local_notify = 0;
        } else {
            $local_notify = 1;
        }
    } elsif ($arg =~ m/^-test/ ) {    
        $local_test = 1;
    } elsif ($arg =~ m/^-e/) {
        setGlobal("utgErrorRate",shift @ARGV);
    } elsif ($arg =~ m/^-g/) {
        setGlobal("utgGenomeSize",shift @ARGV);
    } elsif ($arg =~ m/^-j/) {
        setGlobal("astatLowBound",shift @ARGV);
    } elsif ($arg =~ m/^-k/) {    
        setGlobal("astatHighBound",shift @ARGV);
    } elsif ($arg =~ m/^-(no)?ubs/) {
        setGlobal("utgBubblePopping",1) if ( $1 ne 'no' );
        setGlobal("utgBubblePopping",0) if ( $1 eq 'no' );
    } else {
        $found = 0;
    }

    if ($found) {
        $arg = shift @ARGV;
    }

    return($arg, @ARGV);
}



sub localSetup($) {
    my $numSteps = shift @_;

    asdbInit()  if ( !runningOnGrid());

    #catmap
    if ( -e "$asm.catmap" and !-e "$wrk/$asm.catmap" ) {
      copy("$asm.catmap", "$wrk/$asm.catmap") or die "Could not copy: $asm.catmap\n";
    }

    #seq.features
    if ( -e "$asm.seq.features" and !-e "$wrk/$asm.seq.features" ) {
      copy("$asm.seq.features", "$wrk/$asm.seq.features") or die "Could not copy: $asm.seq.features\n";
    }

    createInvocationScript() && init_prop_file($asm, $numSteps) if (!runningOnGrid());
}



sub localStart ($) {
    my $cmd_name = shift;
    my $props_file = $local_props_file;

    if (! -e "$wrk/log/$cmd_name.started") {
        touch("$wrk/log/$cmd_name.started");

        if (-x $ca_observer ){
            my $exec_cmd = "$ca_observer --appendlog=1 --logfile=$wrk/$ca_log --event=start --name=\"$cmd_name\" --retval=0 --props=$wrk/$props_file -host=`hostname` --message=\"Command with name: '$cmd_name' started\"\n";
            system($exec_cmd);
        }
    }
}



sub localFinish ($) {
    my $cmd_name = shift;
    my $props_file = $local_props_file;

    touch("$wrk/log/$cmd_name.finished");

    if (-x $ca_observer ){
        my $exec_cmd = "$ca_observer --appendlog=1 --logfile=$wrk/$ca_log --event=finish --name=\"$cmd_name\" --retval=0 --props=$wrk/$props_file --host=`hostname` --message=\"Command with name: '$cmd_name' finished\"\n";
        system($exec_cmd);
    }
}



sub localFailure ($) {
    my $msg        = shift @_;
    my $props_file = $local_props_file;

    if (-x $ca_observer ){
        my $exec_cmd = "$ca_observer --appendlog=1 --logfile=$wrk/$ca_log --event=failure --name=\"$0\" --retval=0 --props=$wrk/$props_file --host=`hostname` --message=\"Command with name: '$0' failed\"\n";
        system($exec_cmd);
    }

    open(F, "> $wrk/log/ca.failed");
    print F "$msg\n";
    close(F);
}


sub localPostTerminator($) {
    my $termDir = shift @_;

    link "$termDir/$asm.qc",  "$termDir/$asm.qc.metrics" if (! -e "$termDir/$asm.qc.metrics");
    link "$termDir/$asm.qc.metrics", "$wrk/$asm.qc.metrics";

    my @tempArr = split (/\//, $wrk);
    $request_id = $tempArr[$#tempArr];

    system("cd $wrk && $bin/ametrics.plx -debug 9 $request_id $wrk/$asm.qc.metrics > $wrk/log/ametrics.log 2>&1");
}




# database handle
use IO::File;
use Sys::Hostname;     #Use hostname
use strict; 
use FindBin qw($Bin);

my $dbh = undef;

#Name:   init_db_connection
#Input:  none
#Output: none
#Usage:  This function initializes a connection to the asdb database.
sub init_db_connection() {
    my $sybase = $ENV{SYBASE};
    $ENV{SYBASE} = '/usr/local/packages/sybase'
    	if ( !defined $sybase or !-d $sybase );


    # Try to connect to the database server
    $dbh = DBI->connect(
        "dbi:Sybase:server=SYBTOOLS;packetSize=8092",
        'access', 'access',
        {
            PrintError => 0,
            RaiseError => 0,
            AutoCommit => 1
        }
    );
    $dbh->do("use asdb")
      or die("Failed to open database \'asdb\'");
    $dbh->{InactiveDestroy} = 1;
}



#Name:   asdbInit
#Input:  project name and alias for request id
#Output: none
#Usage:  ASDB initializing for the new request
sub asdbInit () {
    my $project;

    if ( $local_test ) {
    	$local_alias = $asm if (!defined($local_alias);
        $project = "test\n";
    } else {
        $project = $asm
    }


    #First,
    #select service_type_id from Service_Type

    init_db_connection() unless defined $dbh;

    my $query =
"select service_type_id from Service_Type where service_type = 'ASSEMBLY'";


    my $qh = $dbh->prepare($query)
      or die("Cannot prepare $query: " . $dbh->errstr );

    defined $qh->execute()
      or die("Database query \'$query\' failed: " . $dbh->errstr );

    my @row             = $qh->fetchrow();
    my $service_type_id = $row[0];
    $qh->finish();
    $qh = undef;

    #Second
    #insert into Request

    $query =
"insert into Request (status,job_entered,submitter,submitter_host,service_type_id,job_started,job_terminated,submission_type,alias) "
      . "values ('R',getdate(),'"
      . getpwuid($<) . "','"
      . hostname()
      . "',$service_type_id,"
      . 'convert(smalldatetime,getdate()),' . 'NULL,' . '1';

    if ( defined $local_alias ) {
        $query .= ",'$local_alias')";
    }
    else {
        $query .= ",NULL)";
    }

    my $qh = $dbh->prepare($query)
      or die("Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or die("Database query \'$query\' failed: " . $dbh->errstr );

    $qh->finish();
    $qh = undef;

    #Third
    #get generated request id
    $query = "select \@\@identity";


    my $qh = $dbh->prepare($query)
      or die("Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or die("Database query \'$query\' failed: " . $dbh->errstr );

    my @row        = $qh->fetchrow();
    $request_id = $row[0];
    $qh->finish();
    $qh = undef;

    #Fourth
    #update project name in asdb
    $query =
        "insert into Project (request_id, project_name) "
      . "values ($request_id, \"$project\")";

    my $qh = $dbh->prepare($query)
      or die( "Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or die( "Database query \'$query\' failed: " . $dbh->errstr );

    $qh->finish();
    $qh = undef;

    #Fifth
    #Update progress in Process table
    $query =
        "insert into Process (request_id,pid,progress,cluster_id) "
      . "values ($request_id,1,0,0)";

    my $qh = $dbh->prepare($query)
      or die( "Cannot prepare $query: " . $dbh->errstr );
    defined $qh->execute()
      or die( "Database query \'$query\' failed: " . $dbh->errstr );

    $qh->finish();
    $qh = undef;

    print "Your Assembly Console request id is: $request_id\n";
    print "Please check the console for your job's status at\n";
    print "http://aserver.tigr.org:8080/AserverConsole/\n";
    
}

#Name:   createInvocationScript
#Input:  none
#Output: none
#Usage:  Create invocation script
sub createInvocationScript() {
    my $logDir = "$wrk/log";
    my $invocation_script = "$logDir/invocation.info";
    #writing the invocation
    return 0 if ( -e $invocation_script );
        
    my $invo_fh = new IO::File ">$invocation_script"
      or die("Cannot open '$invocation_script'. Error code: $!");
         
    my $cust_name  = getpwuid($<);
    my $hostname   = hostname();
    my $curr_dir = getcwd;
    print $invo_fh "invocation: $invocation\n";
    print $invo_fh "username: $cust_name\n";
    print $invo_fh "hostname: $hostname\n";
    print $invo_fh "userdir: $curr_dir\n";
    close $invo_fh 
      or bail("Error closing '$invocation_script'. Error Code: $!");
    chmod 0755, $invocation_script;

}

#Name:   init_prop_file
#Usage:  Initialize the properties file
sub init_prop_file($$) {
    my $prefix    = shift;
    my $totalcmds = shift;

    my $host = hostname();
    my $pf = $wrk . '/' . $local_props_file;
    #print "Using propsfile: $pf\n";
    my $props_fh = new IO::File ">$pf"
      or die("Could not open props file: '$pf'. Error Code: $!");
    print $props_fh qq~
    request_id=$request_id
    server=SYBTOOLS
    user=access
    password=access
    database=asdb
    host=$host
    clean=9
    prefix=$prefix
    total_commands=$totalcmds
    command_count=0
    ~;
    if ( defined $wrk ) {
        print $props_fh "resdir=$wrk\n";
    }
    if ( $local_notify == 1 && $local_test == 0 ) {
        print $props_fh "email=$local_notify\n";
        print $props_fh "support_email=DLBCIS\n";
    }
    close $props_fh 
      or die("Error closing '$pf'. Error Code: $!");
}

sub copyBack() {

    if ( $local_copy < getGlobal('noCopy') ) {

       my $invocation_script = "$wrk/log/invocation.info";        
       my $invo_fh = new IO::File "<$invocation_script"
	 or die("Cannot open '$invocation_script'. Error code: $!");
       my @invo = <$invo_fh>;

       my $user_dir = undef;
       foreach my $line (@invo) {
       	   chomp($line);
    	   if ( $line =~ /^userdir: (\S+)/ ) {
	     $user_dir = $1;
	   }
       }
       close $invo_fh 
	 or bail("Error closing '$invocation_script'. Error Code: $!");

       if ( defined $user_dir and $user_dir ne '' ) {
	  my $aget = "$Bin/aget.pl";
	  my $agetCmd =  "$aget -mode $local_copy $request_id $user_dir >> $wrk/log/AGET.log";       
	  print "Copying: '$agetCmd'\n";
	  system($agetCmd);
       }       
    }
}

1;
