#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

package libBri;  # is no more.

sub import () {
    print STDERR "This package is obsolete.  Functionality has been moved\n";
    print STDERR "into sim4polish.pm, fasta.pm, and scheduler.pm.  A few\n";
    print STDERR "ESTmapper specific routines moved to ESTmapper/util/\n";
    die;
}

1;
