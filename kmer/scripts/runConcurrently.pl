#!/usr/local/bin/perl

use strict;

#  We use the libBri package of usefull stuff.  It's located in the same place
#  as the ESTmapper.pl script.  That's what FindBin tells us.
#
use FindBin;
use lib "$FindBin::Bin";
use libBri;

&libBri::schedulerSetNumberOfProcesses(4);
&libBri::schedulerSetNumberOfProcesses($ARGV[0]) if (scalar @ARGV > 0);

&libBri::schedulerSetShowCommands(1);
&libBri::schedulerSetShowStatus(1);

while (<STDIN>) {
    &libBri::schedulerSubmit($_);
}

&libBri::schedulerFinish();
