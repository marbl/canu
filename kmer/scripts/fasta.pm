#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

package fasta;

use strict;

$| = 1;

sub import () {
}


######################################################################
#
#  Function to stream a multi-fasta file.  Usage:
#
#  my $header;
#  my $sequence;
#
#  open(F, "< $ARGV[0]");
#  while (!eof(F)) {
#    ($header, $sequence) = nextFastA(*F);
#    print ">$header\n";
#    print "$sequence\n";
#  }
#  close(F);
#
sub nextFastA {
    my ($fh) = @_;
    my $h;
    my $s;

    if (!eof($fh)) {
        my $sep = $/;        #  Save the file line separator

        $/ = "\n";
        $h = <$fh>;

        $h =~ s/^\s+//;
        $h =~ s/\s+$//;

        #chomp $h;

        $/ = ">";
        $s = <$fh>;
        $s =~ s/\s//gm;

        #  Fix things -- the '>' is at the start of the header only if this
        #  is the first sequence read from the file.  Likewise, the sequence
        #  will have the '>' at the end if it is not the last sequence read.
        #
        $h =~ s/^>//g;
        $s =~ s/>$//g;

        $/ = $sep;        #  Restore the file line separator
    }

    return($h, $s);
}
