use strict;
use FindBin qw($Bin);

my $ca_observer = "$Bin/ca_observer.plx";
my $ca_log = "log/ca_observer.log";

sub start ($) {
    if ( !$JCVI ) {
    	return;
    }

    my $cmd_name = shift;
    my $props_file = getGlobal('props_file');

    if ( $JCVI == 1 and -x $ca_observer ){
        touch("$wrk/log/$stepName.started");
        my $exec_cmd = "$ca_observer --appendlog=1 --logfile=$wrk/$ca_log --event=start --name=\"$cmd_name\" --retval=0 --props=$wrk/$props_file -host=`hostname` --message=\"Command with name: '$cmd_name' started\"\n";
        `$exec_cmd`;
    }
}

sub finish ($) {
    if ( !$JCVI ) {
    	return;
    }
    
    my $cmd_name = shift;
    my $props_file = getGlobal('props_file');

    if ( $JCVI == 1  and -x $ca_observer ){
        touch("$wrk/log/$stepName.finished");
        my $exec_cmd = "$ca_observer --appendlog=1 --logfile=$wrk/$ca_log --event=finish --name=\"$cmd_name\" --retval=0 --props=$wrk/$props_file --host=`hostname` --message=\"Command with name: '$cmd_name' finished\"\n";
        `$exec_cmd`;
    }
}

sub failure ($) {
    if ( !$JCVI ) {
    	return;
    }

    my $cmd_name = shift;
    my $props_file = getGlobal('props_file');

    if ( $JCVI == 1  and -x $ca_observer ){
        my $exec_cmd = "$ca_observer --appendlog=1 --logfile=$wrk/$ca_log --event=failure --name=\"$cmd_name\" --retval=0 --props=$wrk/$props_file --host=`hostname` --message=\"Command with name: '$cmd_name' failed\"\n";
        `$exec_cmd`;
	exit;
    }
}


1;
