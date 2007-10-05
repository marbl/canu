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
    if ( !$JCVI ) {
    	return;
    }
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
    if ( !$JCVI ) {
    	return;
    }

    my $project;
    my $alias;

    if ( getGlobal('test') ) {
    	setGlobal('alias',$asm) if ( !defined getGlobal('alias') );
        $project = "test\n";
    } else {
        $project = $asm
    }
    $alias = getGlobal('alias');


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

    if ( defined $alias ) {
        $query .= ",'$alias')";
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
    if ( !$JCVI ) {
    	return;
    }
    
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
    if ( !$JCVI ) {
    	return;
    }
    
    my $prefix    = shift;
    my $totalcmds = shift;

    my $notify    = getGlobal('notify');
    
    my $host = hostname();
    my $pf = $wrk . '/' . getGlobal('props_file');
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
    if ( $notify == 1 && getGlobal('test') == 0 ) {
        print $props_fh "email=$notify\n";
        print $props_fh "support_email=DLBCIS\n";
    }
    close $props_fh 
      or die("Error closing '$pf'. Error Code: $!");
}

sub copyBack() {
    if ( !$JCVI ) {
    	return;
    }
    
    my $copyMode = getGlobal('copy');

    if ( $copyMode < getGlobal('noCopy') ) {

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
	  my $agetCmd =  "$aget -mode $copyMode $request_id $user_dir >> $wrk/log/AGET.log";       
	  print "Copying: '$agetCmd'\n";
	  system($agetCmd);
       }       
    }
}

1;
