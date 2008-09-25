#!/usr/local/bin/perl

use FindBin;
use lib "$FindBin::Bin";
use scheduler;
use strict;

&scheduler::schedulerSetNumberOfProcesses(4);
&scheduler::schedulerSetNumberOfProcesses($ARGV[0]) if (scalar @ARGV > 0);

while (<STDIN>) {
    &scheduler::schedulerSubmit($_);
}

&scheduler::schedulerFinish();
