#!/usr/local/bin/perl

use strict;

my @dateStrings = ( "???", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" );
my %dateStrings;

$dateStrings{"Jan"} = "01";
$dateStrings{"Feb"} = "02";
$dateStrings{"Mar"} = "03";
$dateStrings{"Apr"} = "04";
$dateStrings{"May"} = "05";
$dateStrings{"Jun"} = "06";
$dateStrings{"Jul"} = "07";
$dateStrings{"Aug"} = "08";
$dateStrings{"Sep"} = "09";
$dateStrings{"Oct"} = "10";
$dateStrings{"Nov"} = "11";
$dateStrings{"Dec"} = "12";

if (! -e "logs") {
    system("git log --name-status > logs");
}

#  Update this after each copyright update commit, please.  Best method is to commit
#  the copyright changes -- none of the addCopyrights files, -- update this file and
#  Then commit the addCopyrights files.

my %stoppingCommits;

$stoppingCommits{"6950cb74e302a97673a5ba482b3b8992eea72c37"} = 1;   #  20 AUG 2015 - Initial copyright addition.
$stoppingCommits{"72c27c95d61cb8f37e859c4039456eb2acc5c55b"} = 1;   #  19 NOV 2015 - Second copyright addition.
$stoppingCommits{"b2df5790f77d38cc31fe77a7f65360e02389f92e"} = 1;   #  04 MAR 2016
$stoppingCommits{"1ef335952342ef06ad1651a888f09c312f54dab8"} = 1;   #  18 MAY 2016

open(F, "< logs") or die "Failed to open 'logs': $!\n";

$_ = <F>;  chomp;

my $author;
my $date;

while (!eof(F)) {
    my $commit;

    if (m/^commit\s+(\w+)$/) {
        $commit = $1;
    } else {
        die "Expected commit line, got '$_'\n";
    }

    last  if (exists($stoppingCommits{$commit}));

    $_ = <F>;  chomp;

    if (m/^Merge/) {
        #  Merge commits include an extra line here:  "Merge: a75ed40 717c0b1"
        #  No files change, and we can just process the rest of the entry as normal.
        $_ = <F>;  chomp;
    }

    if      (m/walenz/i) {
        $author = "Brian P. Walenz";
    } elsif (m/koren/i) {
        $author = "Sergey Koren";
    } else {
        print STDERR "Skipping commit from '$_'\n";
        $author = undef;
    }

    $_ = <F>;  chomp;

    if (m/Date:\s+\w+\s+(\w+)\s+(\d+)\s+\d+:\d+:\d+\s+(\d+)/) {
        my $day  = substr("00$2", -2);
        my $mo   = $dateStrings{$1};
        my $year = $3;

        die "Invalid month '$3'\n"  if (! defined($mo));

        $date = "$year$mo$day"

    } else {
        die "Failed to match date in '$date'\n";
    }

    #print STDERR "$commit -- $date -- $author\n";

    $_ = <F>;  chomp;

    while (! m/^commit/) {
        next if (m/^$/);    #  Blank line
        next if (m/^\s+/);  #  Comment line

        if      ($_ =~ m/M\s+(\S+)$/) {
            print "A $1 nihh$date$author\n"   if (defined($author));

        } elsif ($_ =~ m/A\s+(\S+)$/) {
            #  New file, treat as normal.
            print "A $1 nihh$date$author\n"   if (defined($author));

        } elsif ($_ =~ m/D\s+(\S+)$/) {
            #  Deleted file, do nothing,.

        } else {
            print STDERR "$_\n";
        }

    } continue {
        $_ = <F>;  chomp;
    }
}
close(F);
