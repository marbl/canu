#!/usr/local/bin/perl
use strict;
use DBI;
use TIGR::ConfigFile;

my $dbh         = undef;
my $ARUN_CONFIG_FILE = "ARUN.conf";
my $DEFAULT_INSTALL_DIR = "/usr/local/common/ARUN";

my $install_dir = $ENV{'ARUN_INSTALL_DIR'};
$install_dir = $DEFAULT_INSTALL_DIR if ( !defined $install_dir );

my $configFile = $install_dir . '/' . $ARUN_CONFIG_FILE;

print "Config file: '$configFile'\n";
my $arun_cf = new TIGR::ConfigFile($configFile)
    or die("Could not initialize the arun config file object: '$configFile'");

# Default parameters
my $type = $arun_cf->getOption('type','jcvi');
my $server = $arun_cf->getOption('server','jcvi');
my $user = $arun_cf->getOption('db_user','jcvi');
my $pass = $arun_cf->getOption('db_pass','jcvi');
my $database = $arun_cf->getOption('database','jcvi');

#This function initializes a connection to the asdb database.
sub init_db_connection() {
	# Try to connect to the database server
	$dbh = DBI->connect(
		"dbi:Sybase:server=$server;packetSize=8092",
		"$user", "$pass",
		{
			PrintError => 0,
			RaiseError => 0,
			AutoCommit => 1
		}
	);
	$dbh->do("use $database")
	  or die("Failed to open database \'$database\'");

	$dbh->{InactiveDestroy} = 1;
}
MAIN:
{
   my $request_id  = $ARGV[0];

   my $start_query =
      'update Request set '.
      "status = \'R\' where request_id = $request_id";
   print("Executing the status query $start_query",3);
   if ( !defined $dbh ) {
	   init_db_connection();
   }
   my $qh = $dbh->prepare($start_query) or
   die("Cannot prepare $start_query: " . $dbh->errstr);
   if (! defined $qh->execute()) {
      die("Database query \'$start_query\' failed: " .
    			  $dbh->errstr);
   }
   $dbh->disconnect();
   $dbh = undef;
}

