#!/usr/local/bin/perl

use strict;
use DBI;
use IO::File;
use TIGR::Foundation;
use Sys::Hostname;
use File::Copy;
use File::Find;
use HTTP::Headers;
use HTTP::Request;
use URI::URL;
use LWP::UserAgent;
use MIME::Lite;

my $HELPTEXT = qq~
 
 This is the observer script that logs the progress of an Assembly job 
 in the database. It is called directly by the Workflow system

~;

my $REQUEST_URL = "http://HOST:18080/request/req_table?request_id=REQUEST_ID";

my $CMD_LOG = "command.out";
my $tf_object = new TIGR::Foundation();
my $name = undef;
my $event = undef;
my $inst_filename = undef;
my $prop_filename = undef;
my $message = undef;
my $dbh = undef;
my $request_id = undef;
my $server = undef;
my $user = undef;
my $password = undef;
my $database = undef;
my $host = undef;
my $tot_count = undef;
my $cmd_count = undef;
my $resdir = undef;
my $marshalling_base = undef;
my $outdir = undef;
my $email = undef;
my $support_email = undef;
my $text = undef;
my $prefix = undef;

our @DEPEND = ("TIGR::Foundation");
our $REVISION = (qw$Revision: 1.1 $)[-1];
our $VERSION = '2.2';
our $VERSION_STRING = "$VERSION (Build $REVISION)";

#This function initializes the connection to the database
sub init_db_connection() {
 
   print("Establishing a connection to the database");
   my $type = "Sybase";
   # Try to connect to the database server
   $dbh = DBI->connect("dbi:$type:server=$server;packetSize=8092",
                         $user, $password,
                         { PrintError => 0,
                           RaiseError => 0,
                           AutoCommit => 1  
                         }
                      );
   $dbh->do("use $database") or 
      die("Failed to open database \'$database\'");
}

#Name:   updateErrorStatus($keyword)
#Input:  the keyword indicating the message to be put in the Messages table
#Usage:  This function updates the Request table and puts the appropriate
#        message in the Messages table to indicate that there was an error 
#        while processing the request.
sub UpdateErrorStatus($;$) {
   my $keyword = shift;
   my $params = shift;
   
   if(!defined $dbh) {
      init_db_connection();
   }
   if(defined $dbh) {
      my $query = "select message_code from Message_Template where template like ".
                  "\'%$keyword%\'"; 
   
      print("Running query $query ...");
      my $qh = $dbh->prepare($query) or die("Cannot prepare $query: " . 
                                          $dbh->errstr);
      if(! defined $qh->execute())
      {
         die("Database query \'$query\' failed: " . $dbh->errstr);
      }
  
      my $code = $qh->fetchrow();
      my $status = "E";
   
      $qh = undef;
      $query = "update Request set job_terminated = convert(smalldatetime,getdate()), ".
            "status = \'$status\' where request_id = $request_id";   
   
      $qh = $dbh->prepare($query) or die("Cannot prepare $query: " .$dbh->errstr);
      if(! defined $qh->execute())
      {
         die("Database query \'$query\' failed: " . $dbh->errstr);
      }
      $qh = undef;

      if(defined $params) {
         $query = "exec putMessage \@message_code=$code, \@request_id=$request_id, ".
                  "\@params=\"$params\"";
      }
      else {
	 $query = "exec putMessage \@message_code=$code, \@request_id=$request_id";
      }

      $qh = $dbh->prepare($query) or die("Cannot prepare $query: " .$dbh->errstr);
      if(! defined $qh->execute())
      {
         die("Database query \'$query\' failed: " . $dbh->errstr);
      } 
      $qh = undef;
   }
   $dbh->disconnect();
   $dbh = undef;
 } 


#This method updates the Progress table in the asdb database 
#with the percentage of assembly completion
sub putProgress($) {
   my $pf = shift;

   my $new_count = $cmd_count + 1;
   my $new_progress = ($new_count/$tot_count)*100;
   $new_progress = sprintf("%2.0f", $new_progress);
   my $pid = getppid;
   my $process_query = "update Process set progress=$new_progress ".
                       "where request_id = $request_id";
     
   print("Executing the process query $process_query");
   
   my $qh = $dbh->prepare($process_query) or 
   die("Cannot prepare $process_query: " . $dbh->errstr);
   
   if (! defined $qh->execute()) {
      die("Database query \'$process_query\' failed: " . $dbh->errstr);
   }
   $pf->print("command_count=$new_count\n");
   if($tot_count == $new_count) {

      my $status = undef;
      #Updating the Request table for the finished job
      $status = "F";
      my $status_query = 
        "update Request set job_terminated = convert(smalldatetime,getdate()), ".
        "job_retrieved = convert(smalldatetime,getdate()), ".
        "status = \'$status\' where request_id = $request_id";
   
      print("Executing the status query $status_query");
   
      my $qh = $dbh->prepare($status_query) or 
      die("Cannot prepare $status_query: " . $dbh->errstr);
   
      if (! defined $qh->execute()) {
         die("Database query \'$status_query\' failed: " . 
                           $dbh->errstr);

      }
      $dbh->disconnect();
      $dbh = undef;
      
   } 
}

sub sendMail() {
   if(!defined $dbh) {
      init_db_connection();
   }
   my $stats_query = 
      "select status, submitter from  Request where ".
      "request_id = $request_id";
   
   print("Executing the query $stats_query");
   
   my $qh = $dbh->prepare($stats_query) or 
      die("Cannot prepare $stats_query: " . $dbh->errstr);
   
   if (! defined $qh->execute()) {
      die("Database query \'$stats_query\' failed: " . $dbh->errstr);
   }
   my @row = $qh->fetchrow();
   my $status = $row[0];
   my $submitter = $row[1];
   
   $qh->finish;
   $qh = undef;
   
   $dbh->disconnect();
   $dbh = undef;

   my $subject = undef;
   
   if($status eq "F") {
       $subject = "Assembly request $request_id completed";
   }
   if($status eq "E") {
       $subject = "Assembly request $request_id had errors";
   }
   
   my $body = undef;
   
   if(!$text) {
      my $request_url = $REQUEST_URL;
      $request_url =~ s/HOST/$host/;
      $request_url =~ s/REQUEST_ID/$request_id/;
      my $url = new URI::URL($request_url);
      my $req = new HTTP::Request("GET", $url);
      my $ua = new LWP::UserAgent;
      my $resp = $ua->request($req);
      if ($resp->is_success) {
         $body = $resp->content."\n";
      }
      else {
         my $message = $resp->message;
         print("Could not send mail from $REQUEST_URL : $message"); 
         $text = 1;
      }   
   }

   if($text) {
      if($status eq "E") {
         $body = "Your request with id $request_id failed.\n".
              "Please see the error files at $marshalling_base/$request_id \n".
              "for error tracking. Also check for error messages from the \n".
              "job by visiting the AserverConsole at /http://aserver.tigr.org:18080/aserver.\n";
      }  
      if($status eq "F") {
         $body = "Your request with id $request_id was successful.\n";
         if(-e $resdir) {
	    $body .= "The assembly results can be seen at $resdir.\n";
	 }
      }
   }

   if($email) {
      if($text) {
         use Mail::Mailer;
         my $type = 'sendmail';
         my $mailprog = Mail::Mailer->new($type);
         # mail headers to use in the message
         my %headers = (
            'To' => "$submitter\@jcvi.org",
            'From' => 'aserver@jcvi.org',
            'Subject' => $subject 
         );
   
         $mailprog->open(\%headers);
         print $mailprog $body;
         $mailprog->close;
      }
      else {
         my $msg = MIME::Lite->new(
                        From    =>"aserver\@jcvi.org",
                        To      =>"$submitter\@jcvi.org",
                        Subject =>"$subject",
                        Type    =>"multipart/related"
                       );
         

         $msg->attach(Type => "text/html",
                      Data => $body);
         $msg->send();
      }
   }

   if(($status eq "E") && (defined $support_email)) {
      my $subject = "Assembly request $request_id had errors";;
      my $body = undef;
      
      $body = "The request with id $request_id failed.\n".
         "Please see the error files at $marshalling_base/$request_id \n".
         "for error tracking. Also check for error messages from the \n".
         "job by visiting the Aserver Console at /http://aserver.tigr.org:18080/aserver.\n";
   
      use Mail::Mailer;
      my $type = 'sendmail';
      my $mailprog = Mail::Mailer->new($type);
      # mail headers to use in the message
      my %headers = (
         'To' => "$support_email\@jcvi.org",
         'From' => "$submitter\@jcvi.org",
         'Subject' => $subject 
      );
   
      $mailprog->open(\%headers);
      print $mailprog $body;
      $mailprog->close;
   } 
}      

sub bailOut($) {
   my $message = shift;
   #ending email on error
   if($email) {
      print("Sending mail to the job submitter ".
                     "when job errors out");
      sendMail();
   }
   die("$message");
}

MAIN: 
{
   my $name = undef;
   my $id = undef;
   my $time = undef;
   my $event = undef;
   my $inst_filename = undef;
   my $prop_filename = undef; 
   my $fail_message = undef;
   my $hostname = undef;
   my $message = undef;
   my $retval = undef;
   my $pf = undef;

   $tf_object = new TIGR::Foundation;
   # add dependency information
   $tf_object->addDependInfo(@DEPEND);
   # add version information
   $tf_object->setVersionInfo($VERSION_STRING);
   # add help information
   $tf_object->setHelpInfo($HELPTEXT);

   #Read command line options
   my $result = $tf_object->TIGR_GetOptions (
                                    'name=s',   \$name,
                                    'ID=i',     \$id,
                                    'time=s',   \$time,
       				    'event=s',  \$event,
       				    'file=s',   \$inst_filename,
                                    'props=s',  \$prop_filename,
                                    'message=s', \$message, 
                                    'retval=i',  \$retval,
                                    'fail=s',    \$fail_message,
                                    'host=s',   \$hostname
				  );
   
   if ( ! defined ( $result ) ) {
      die("The options could not be read");
   }
   #Set the logfile for the observer. This file is appended
   $tf_object->logAppend(1);
   $tf_object->setDebugLevel(1);
   my $prog_name = $tf_object->getProgramInfo("name");
   my $invocation = $tf_object->getProgramInfo("invocation");
   print("START: $prog_name $invocation");
   
   #read the properties file to get important information about the job
   if((defined $prop_filename) && (-e $prop_filename)) {
      $pf = new IO::File("$prop_filename") or 
              die("Failed to open  $prop_filename($!)");
      my $line = undef;
      my @props_arr = ();

      while(defined ($line = <$pf>)) {
         chomp($line);
         @props_arr = split("=", $line);
           
         if($props_arr[0] =~ /request_id/) {
	    $request_id = $props_arr[1];
            print("The request id is $request_id");
	 }
         
         if($props_arr[0] =~ /server/) {
	    $server = $props_arr[1];
            print("The server is $server");
	 }
  
         if($props_arr[0] =~ /user/) {
	    $user = $props_arr[1];
            print("The user is $user");
	 }
    
         if($props_arr[0] =~ /password/) {
	    $password = $props_arr[1];
            print("The password is $password");
	 }
        
         if($props_arr[0] =~ /database/) {
	    $database = $props_arr[1];
            print("The database is $database");
	 }
 
         if($props_arr[0] =~ /host/) {
	    $host = $props_arr[1];
            print("The host is $host");
	 }

         if($props_arr[0] =~ /resdir/) {
	    $resdir = $props_arr[1];
            print("The resdir is $resdir ");
	 }
  
         if($props_arr[0] =~ /marshalling_base/) {
	    $marshalling_base = $props_arr[1];
            print("The marshalling_base is $marshalling_base");
	 }
         
         if((defined $marshalling_base) && (defined $request_id)) {
            $outdir = "$marshalling_base/$request_id";
            print("the outdir is $outdir");
	 }

         if($props_arr[0] =~ /email/) {
	    $email = $props_arr[1];
            print("The email is $email");
	 }

         if($props_arr[0] =~ /support/) {
	    $support_email = $props_arr[1];
            print("The support email is $support_email");
	 }

         if($props_arr[0] =~ /text/) {
	    $text = $props_arr[1];
            print("The text is $text");
	 }
 
         if($props_arr[0] =~ /prefix/) {
	    $prefix = $props_arr[1];
            print("The prefix is $prefix");
	 }

         if($props_arr[0] =~ /total_commands/) {
	    $tot_count = $props_arr[1];
            print("The total number of commands is $tot_count");
	 }
  
         if($props_arr[0] =~ /command_count/) {
	    $cmd_count = $props_arr[1];
         }

         shift @props_arr;
      }
   }
   if(defined $cmd_count) {
      print("The command count is $cmd_count");
   }        
   close $pf;
   
   #Initialize a connection to the database
   init_db_connection();
 
   if(($event eq "start") && (defined $hostname)) { 
                                                   #if a program started log the 
                                                   #step in the 
                                                   #AdminLog table in asdb
      chomp $name;
      my $status = "$name started";
      my $user = getpwuid($<);
      if (defined($hostname)) {
	  $status .= " on $hostname";
      }
             
      if(!defined $user) {
         die("The user is not defined");
      }
      my $status_query = "exec $database..putAdminLog '$status','$user',$request_id";
      print("Running query $status_query ...");
      my $qh = $dbh->prepare($status_query) or 
         die("Cannot prepare $status_query: " .$dbh->errstr);
      
      if(! defined $qh->execute())
      {
         die("Database query \'$status_query\' failed: " . $dbh->errstr);
      }

      $qh = undef;
      
      my $request_query = "update Request set job_host = \"$hostname\" ".
               "where request_id = $request_id";
      print("Running query $request_query ...");
      $qh = $dbh->prepare($request_query) or 
         die("Cannot prepare $request_query: " .$dbh->errstr);
      
      if(! defined $qh->execute())
      {
         die("Database query \'$request_query\' failed: " . 
                              $dbh->errstr);
      }
      $qh = undef;
   }
  
   if(($event eq "failure") || ((defined $retval) && ($retval != 0))){ #if the command failed update the Request 
                                                #and Messages table in asdb 
      UpdateErrorStatus("CA pipeline", $name);
      bailOut("ERROR: The $name program from the CA pipeline was not successful");
   }
   
   if($event eq "finish") { #if the command finished Update the Process table 
                            #and the RequestTable
      $pf = new IO::File(">> $prop_filename") or 
              die("Failed to open  $prop_filename($!)");
      putProgress($pf);
      close($pf);
   }

   if(defined $dbh) {
      $dbh->disconnect();
   }

   exit 0;   
}
