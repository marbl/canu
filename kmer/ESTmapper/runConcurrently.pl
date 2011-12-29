#!/usr/local/bin/perl

use strict;

use FindBin;
use lib "$FindBin::Bin/../lib";
use scheduler;

&scheduler::schedulerSetNumberOfProcesses(4);
&scheduler::schedulerSetNumberOfProcesses($ARGV[0]) if (scalar @ARGV > 0);

while (<STDIN>) {
    &scheduler::schedulerSubmit($_);
}

&scheduler::schedulerFinish();
