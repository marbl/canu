#!/usr/bin/perl

#  A basic (for now) script to diff to qc files.  In the near future,
#  this will be extended to decide if things changed significantly or
#  not.

use strict;

#  Set to 1 if there are differences.
my $diffs = 0;

if (scalar(@ARGV == 0)) {
    die "usage: $0 new.qc old.qc\n";
}
if (! -e $ARGV[0]) {
    die "Can't find new.qc '$ARGV[0]'\n";
}
if (! -e $ARGV[1]) {
    die "Can't find old.qc '$ARGV[1]'\n";
}

system("ls -l $ARGV[0]");
system("ls -l $ARGV[1]");
print "\n";

system("diff --side-by-side -W 84 -b $ARGV[0] $ARGV[1]");
$diffs = $? >> 8;

exit($diffs);
