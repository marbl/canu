# This file is part of A2Amapper.
# Copyright (c) 2005 J. Craig Venter Institute
# Author: Brian Walenz
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

use strict;
use Config;  #  for @signame

#  Utility to run a command and check the exit status
#
sub runCommand {
    my $cmd = shift @_;

    print STDERR "$cmd\n";

    my $rc = 0xffff & system($cmd);

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!
    #
    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    }

    my $error = "ERROR: $cmd\n        failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } elsif ($rc > 0x80) {
        $rc >>= 8;
        $error .= "exit status $rc\n";
    } else {
        if ($rc & 0x80) {
            $rc &= ~0x80;
            $error .= "coredump from ";
        }
        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc]\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}

1;
