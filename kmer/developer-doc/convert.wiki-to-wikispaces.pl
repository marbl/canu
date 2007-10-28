#!/usr/bin/perl

use strict;

open(F, "< estmapper.wiki")       or die "Failed to open 'estmapper.wiki'\n";
open(O, "> estmapper.wikispaces") or die "Failed to open 'estmapper.wikispaces'\n";

my $tableline;

while (<F>) {
    chomp;

    s!<i>!//!g;
    s!</i>!//!g;

    s!<b>!**!g;
    s!</b>!**!g;

    s!<tt>!{{!g;
    s!</tt>!}}!g;

    s!<br>!!g;

    s!<pre>!{{!g;
    s!</pre>!}}!g;

    s!<sub>!_!g;
    s!</sub>!!g;

    #  Special cases for the big table of options
    next if (m!^\{\|!);
    next if (m!^\|\}!);
    next if (m!^\|-!);

    if (m!^\|valign="top"\|(.*)$!) {
        $tableline = "|| $1";
        next;
    }

    if (defined($tableline)) {
        if (m!^\|(.*)!) {
            $_ = "$tableline || $1 ||";
        } else {
            die "bad table\n";
        }
        undef $tableline;
    }

    if      (m/^;(.*)$/) {
        print O "|| $1 || ";
    } elsif (m/^:(.*)$/) {
        print O "$1 ||\n";
    } elsif (m/=(=+)\s+(.*)\s+(=+)=$/) {
        print O "$1 $2 $3\n";
    } else {
        print O "$_\n";
    }
}

close(F);
close(O);

print STDERR "Edit estmapper.wikispaces to fix the big option table, search for '*{' and merge those lines.\n";
