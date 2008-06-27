#!/usr/local/bin/perl -w

use strict;
use TIGR::Foundation;
use DBI;
use TIGR::ConfigFile;
use IO::File;
use TIGR::SP;

my $HELPTEXT = qq~
ametrics.plx - a program that reads the .metrics file from the CA pipeline
and uploads metrics information into the asdb database.

ametrics.plx [options] <request_id> <metricsfile>
   request_id  - The request_id of an Assembly request that was
                 received by Aserver.
   metricsfile - The file that is in INI format and contains the
                 metrics information.

Options:
   -__conf__ <aserver_config>  Specify a config file for the aserver(by
                        default use the config file in
                        /usr/local/packages/aserver/etc/aserver.conf).

Exit Codes:
   0   The ametrics program completed successfully;
   1   There were problems while executing ametrics.
~;

my $tf= new TIGR::Foundation;

my $VERSION = " Version 2.00 (Build " . (qw/$Revision: 1.2 $/ )[1] . ")";
my @DEPENDS =
(
   "TIGR::Foundation",
   "TIGR::SP",
   "TIGR::ConfigFile"
);

$tf->addDependInfo(@DEPENDS);
$tf->setHelpInfo($HELPTEXT);
$tf->setVersionInfo($VERSION);

my $METRICS_MAX = 4000;
my $REQUEST = "Request";
my $METRICS  = "Metrics";
my $METRICS2 = "Metrics2"; # For metric values greater than $METRICS_MAX

my $METRIC_TAGS = "Metric_Tags";
my $METRIC_SECTIONS = "Metric_Sections";
my $METRIC_TEMPLATE = "Metric_Template";
my $dbh = undef;
my $BCP_COUNT = 10000;
my $dbserver = undef;   # database server
my $dbtype = undef;      # database type
my $username = undef;
my $password = undef;
# The executable for bcp
my $BCP_EXEC = "/usr/local/packages/sybase/OCS/bin/bcp";
my $CLEAN_METRICS_SP = "clean_metrics_sp";
# default config file for the aserver
my $aserver_config_file = "/usr/local/packages/aserver/etc/aserver.conf";

# hashes that contain metrics information
my %inserted_sections = ();
my %inserted_options = ();
my %inserted_metrics = ();

# cleanTable - This function cleans any records in the Metrics and Metrics2
# tables with the current request_id. The function returns 1 on success.
sub cleanTable($) {
   my $request_id = shift;
   my $sp_obj = new TIGR::SP($dbh) ||
      $tf->bail("No TIGR::SP object created");

   my $sp_name = $CLEAN_METRICS_SP;
   my $sp_params = undef;
   my $rc = undef;
   $tf->logLocal("Cleaning the Metrics tables", 1);

   $sp_params = $request_id;
   $rc = $sp_obj->execute($sp_name, $sp_params);

   my $error = undef;
   if (!defined($rc) ) {# Check fatal error
      # check errors
      my @errs = $sp_obj->error();
      my $errors = join("\n", @errs);
      $tf->bail("Execution of procedure \'$sp_name\' failed because of ".
              "the following errors: $errors");
   }
   # get status
   $error = $sp_obj->status();

   if(defined $error) {
      # collect errors due to stored procedure execution
      my @errs = $sp_obj->error();
      my $errors = join("\n", @errs);

      if ($error < 0) {
         $tf->logError("A Sybase related error occurred");
         $tf->bail("The errors are $errors");
      }
      elsif ($error != 0) {
         $tf->logError("Stored procedure $sp_name failed");
         $tf->bail("The errors are $errors");
      }
   }
   $sp_obj->finish();
   return 1;
}

# This function gets the section, option and metrics information
# stored in the Metric_Sections, Metric_Tags and Metric_Template
# tables
sub getInsertedInfo() {
   my $query = "select s.section_id, s.section, t.tag_id, t.tag, ".
               "p.metric_id from $METRIC_SECTIONS s, $METRIC_TAGS t, ".
               "$METRIC_TEMPLATE p where s.section_id = p.section_id ".
               "and t.tag_id = p.tag_id";

   $tf->logLocal("Running the query $query", 3);
   my $qh = $dbh->prepare($query) or
      $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }

   while (my @row = $qh->fetchrow()) {
      my $section_id = $row[0];
      my $section = $row[1];
      my $tag_id = $row[2];
      my $tag = $row[3];
      my $metric_id = $row[4];

      if((defined $section) && ($section ne "NULL") &&
         ($section ne "")) {
	  $inserted_sections{$section} = $section_id;
      }

      if((defined $tag) && ($tag ne "NULL") &&
         ($tag ne "")) {
	  $inserted_options{$tag} = $tag_id;
      }

      if(defined $metric_id) {
         $inserted_metrics{$section_id}->{$tag_id} = $metric_id;
      }
   }

   $qh->finish();
   $qh= undef;
}

# This function inserts a new section in the Metric_Sections table. It takes
# a new section name and inserts it in the Metric_Sections table returning the
# section_id for the new section.
sub insertSection($) {
   my $section = shift ;
   my $section_id = undef;
   # insert query
   my $query = "insert into $METRIC_SECTIONS (section) values (\"$section\")";
   $tf->logLocal("Running the query $query", 3);
   my $qh = $dbh->prepare($query) or
     $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }
   $qh->finish();
   $qh= undef;

   # get the new section_id
   $query = "select \@\@identity";

   $tf->logLocal("Running the query $query", 3);
   $qh = $dbh->prepare($query) or
            $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }
   my @row = $qh->fetchrow();
   $section_id = $row[0];
   if($section_id == 0) {
      $section_id = undef;
   }
   return $section_id;
}

# This function inserts a new option in the Metric_Tags table. It takes
# a new option name and inserts it in the Metric_Tags table returning the
# option_id for the new section.
sub insertOption($) {
   my $option = shift ;
   my $option_id = undef;
   #insert query
   my $query = "insert into $METRIC_TAGS (tag) values (\"$option\")";

   $tf->logLocal("Running the query $query", 3);
   my $qh = $dbh->prepare($query) or
     $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }
   $qh->finish();
   $qh= undef;

   # get the new option_id
   $query = "select \@\@identity";

   $tf->logLocal("Running the query $query", 3);
   $qh = $dbh->prepare($query) or
            $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }
   my @row = $qh->fetchrow();
   $option_id = $row[0];
   if($option_id == 0) {
      $option_id = undef;
   }
   return $option_id;
}

# This function takes in the value of a metrics and returns the datatype
# to which it belongs. Valid datatypes are String, Date, Integer, Double,
# Float, Other, StringList, DateList, IntegerList, DoubleList, FloatList and
# OtherList.
sub getDataType($) {
   my $value = shift;
   # get the value for the option
   my @values = split(",",$value);

   my $data_type = "Other";
   if (defined($values[0])) {
      $data_type = "String" if ($values[0] =~ /\D/); # nondigits

      if (($data_type eq "String") &&
	  ($values[0] =~ /(\d\d):(\d\d):(\d\d)/)) { # hh:mm:ss
         my ($hh, $mm, $ss) = split(":",$values[0]);
	 if (($hh<24) && ($mm<60) && ($ss<60)) {
	    $data_type = "Date";
	 }
      }
      else {
         @values = split(" ",join(" ",@values));

         foreach my $val (@values) {
	    if ($val =~ /^-?\d+$/) { # an integer
	       $data_type = "Integer";
	       if ($val =~ /^\d{5,}$/) { # whole big number
                  $data_type = "Double";
		  last;
	       }
	    }
            elsif ($val =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)$/) {
               $data_type = "Float";
	       last;
	    }
            else {
               $data_type = "String";
	       last;
	    }
	 }
      }
      $data_type .= "List" if (scalar(@values) > 1);
   }
   return $data_type;
}

# This function inserts a new record in the Metric_Template table. It takes
# in the section_id, option_id, data_type and service_type_id for the new record.
# It returns the metric_id for the new record inserted.
sub insertMetricsTemplate($$$$){
   my $section_id = shift;
   my $option_id = shift;
   my $data_type = shift;
   my $service_type_id = shift;
   #insert new record
   my $query = "insert into $METRIC_TEMPLATE (section_id,tag_id,data_type,".
               "service_type_id) values ($section_id,$option_id,".
               "\"$data_type\",$service_type_id)";

   $tf->logLocal("Running the query $query", 3);
   my $qh = $dbh->prepare($query) or
     $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }
   $qh->finish();
   $qh= undef;

   # get the new metric_id
   $query = "select \@\@identity";

   $tf->logLocal("Running the query $query", 3);
   $qh = $dbh->prepare($query) or
            $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }
   my @row = $qh->fetchrow();
   my $metric_id = $row[0];
   return $metric_id;
}

#This function bcps metrics information into the Metrics and Metrics2 tables.
#It takes in reference to a hash containing metrics information, a request_id,
#the table name for the bcp and the bcp file name. The data is inserted in
#batches of a 10000 records.
sub bcpRecords($$$$) {
   my $bcp_metrics_ref = shift;
   my %bcp_metrics = %$bcp_metrics_ref;
   my $request_id = shift;
   my $bcp_table = shift;
   my $bcp_file = shift;
   my $bcp_output_file = "bcp_out";
   my $addbcp = new IO::File(">>$bcp_file")  or
                 $tf->bail("Cannot open $bcp_file ($!)");

   my $bcp_file_count = 0;
   my $bcp_count_file = "$bcp_file$bcp_file_count";
   my $bcp = new IO::File(">$bcp_count_file")  or
                 $tf->bail("Cannot open $bcp_count_file ($!)");

   my $record = undef;
   my $record_count = 0;

   my @metrics_arr = keys(%bcp_metrics);
   my $metric_id = undef;

   foreach $metric_id (@metrics_arr) {
      $record_count++;
      my $metric_val = $bcp_metrics{$metric_id};
      $bcp->print("$metric_id\t$request_id\t$metric_val\n");
      $addbcp->print("$metric_id\t$request_id\t$metric_val\n");

      if($record_count == $BCP_COUNT ) {
         close($bcp);
         my $bcp_command = "$BCP_EXEC $bcp_table in ".
                "$bcp_count_file ".
                "-c -U $username -P $password -S $dbserver > $bcp_output_file";
         $tf->logLocal("Running \'$bcp_command\' ...", 3);

         my $bad = $tf->runCommand($bcp_command);
         $tf->bail("Command failed: \'$bcp_command\' ($!)") if ($bad);

         $record_count = 0;
         $bcp_file_count++;
         $bcp_count_file = "$bcp_file$bcp_file_count";
         $bcp = new IO::File(">$bcp_count_file")  or
                   $tf->bail("Cannot open $bcp_count_file ($!)");
      }
   }
   close $bcp;
   close($addbcp);

   my $bcp_command = "$BCP_EXEC $bcp_table in ".
          "$bcp_count_file ".
          "-c -U $username -P $password -S $dbserver > $bcp_output_file";
   $tf->logLocal("Running \'$bcp_command\' ...", 3);

   my $bad = $tf->runCommand($bcp_command);
   $tf->bail("Command failed: \'$bcp_command\' ($!)") if ($bad);
   my $debug = $tf->getDebugLevel();
   for(my $i=0; $i <= $bcp_file_count; $i++) {
      unlink "$bcp_file$i";
   }
   unlink $bcp_output_file;

   unlink $bcp_file unless((defined $debug) && ($debug > 0));
}

MAIN:
{
   my %metricVals = ();
   my %big_metricVals = ();
   my $metric_id = undef;

   $tf->TIGR_GetOptions('__conf__=s',  \$aserver_config_file);

   #creating the config file object
   my $acf = new TIGR::ConfigFile($aserver_config_file);

   if(!defined $acf) {
      $tf->bail("Could not initialize config file object");
   }
   $dbtype = "Sybase";
   $dbserver = $acf->getOption('server');
   $username = $acf->getOption('db_user');
   $password = $acf->getOption('db_pass');
   my $db = $acf->getOption('database');

   my $sybase = $ENV{SYBASE};
   $ENV{SYBASE} = '/usr/local/packages/sybase'
       if ( !defined $sybase or !-d $sybase );

   # now we try logging into the specified database
   $dbh = DBI->connect("dbi:$dbtype:server=$dbserver;packetSize=8092",
                        $username, $password,
                        { PrintError => 0,
                          AutoCommit => 1
                        });
   if (! defined $dbh) {
      $tf->bail("ERROR: Connection to server '$dbserver' failed: " .
                  $DBI::errstr);
   }

   $dbh->do("use $db") or
      $tf->bail("ERROR: Failed to open database '$db': " . $dbh->errstr());

   my $request_id = $ARGV[0];
   my $metrics_file = $ARGV[1];

   my $query = "select service_type_id from $REQUEST where ".
               "request_id = $request_id";

   $tf->logLocal("Running the query $query", 3);
   my $qh = $dbh->prepare($query) or
     $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }
   my @row = $qh->fetchrow();
   my $service_type_id = $row[0];

   $qh->finish();
   $qh= undef;

   #creating the metrics file object
   my $cf = new TIGR::ConfigFile($metrics_file);

   if(!defined $cf) {
      $tf->bail("Could not initialize metrics file object");
   }

   #cleaning the existing metrics for this request_id
   if(cleanTable($request_id) == 1) {
      $tf->logLocal("The Metrics tables have been ".
              "successfully cleaned", 1);
   }

   # get all the sections in the metrics file
   my $sections = $cf->getSections();

   $tf->logLocal("Get inserted information from the $METRIC_SECTIONS, ".
           "$METRIC_TAGS and $METRIC_TEMPLATE tables", 1);
   getInsertedInfo();

   my $section = undef;

   foreach $section (@$sections) {
      my $section_id = undef;
      if(!defined($inserted_sections{$section})) { #if the section is not
                                                   #present insert it in
                                                   #Metric_Sections
         $section_id = insertSection($section);
      }
      else {
         $section_id = $inserted_sections{$section};
      }

      my $options_ref = $cf->getOptionNames($section);
      my @options = @$options_ref;
      my $option = undef;

      foreach $option (@options) { #if the option is not
                                   #present insert it in Metric_Tags
         my $option_id = undef;

         if(!defined($inserted_options{$option})) {
            $option_id = insertOption($option);
         }
         else {
	    $option_id = $inserted_options{$option};
         }

         my $value = $cf->getOption($option, $section);

         if((!defined($inserted_sections{$section})) ||
            (!defined($inserted_options{$option}))) { #if either the section or
                                  #option is new insert it in Metric_Template

            # Get the data type for the metric value
            my $data_type = getDataType($value);
            $tf->bail("The datatype is not defined for the value in ".
                    "option $option of section $section")
	       if(!defined $data_type);

            # insert a new record in the Metric_Template table
            $metric_id = insertMetricsTemplate($section_id, $option_id, $data_type,
						  $service_type_id);
         }
         else {
	    $metric_id = $inserted_metrics{$section_id}->{$option_id};
         }

         if($value =~ /^\s*$/) {
	    $value = " ";
         }
         # insert big metric values in the Metrics2 table and the small values
         # in the Metrics tables.
         if((length($value)) > $METRICS_MAX) {
            $big_metricVals{$metric_id} = $value;
	 }
	 else {
            $metricVals{$metric_id} = $value;
	 }
      }
   }

   my @big_metrics = keys(%big_metricVals);
   my @metrics = keys(%metricVals);
   my $bcp_out = undef;

   # bcp the records in the Metrics tables.
   if((scalar(@big_metrics)) > 0) {
      $bcp_out = "big_bcp_file";
      bcpRecords(\%big_metricVals, $request_id, "asdb..Metrics2", $bcp_out);
   }

   if((scalar(@metrics)) > 0) {
      $bcp_out = "bcp_file";
      bcpRecords(\%metricVals, $request_id, "asdb..Metrics", $bcp_out);
   }

   $query = "update $REQUEST set metrics=1 where request_id = $request_id";

   $tf->logLocal("Running the query $query", 3);
   $qh = $dbh->prepare($query) or
      $tf->bail("Cannot prepare $query: " . $dbh->errstr);

   if (! defined $qh->execute()) {
      $tf->logError("Database query \'$query\' failed: " . $dbh->errstr);
      my $ef = (defined $tf->getErrorFile())? $tf->getErrorFile() : "";
      $tf->bail("Error: Failed date query from Sybase, see $ef ".
              "for details. Exiting...");
   }

   $qh->finish();
   $qh= undef;

}

