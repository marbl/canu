#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;
use Config;  #  for @signame
use Time::Local;

#  fetch
#     update the local repositories from github.
#
#  checkout date
#     Date must be of the form yyyy-mm-dd-hhmm.  This will checkout
#     the assembler current at that time.  If not present, the current
#     time is assumed.  A change log is generated from the last
#     date present.
#
#  build date
#     Builds the assembler checked out into 'date'.
#
#  assemble date specFile .... email ...
#     Launches canu on each specFile, running it in directory 'date'.
#     Each specFile must have a unique name; the assembly is named
#     after the specFile.
#     At the end, a diff is generated to POINTERS/$prefix.last and
#     POINTERS/$prefix.reference.  If the assembly finished, POINTERS/$prefix.last
#     is updated.
#
#  submit date
#     Submit ourself for execution at date.
#

my $oper    = shift @ARGV;

my $site    = undef;
my $wrkdir  = undef;
my $gitrepo = undef;

my $tz      = `date +%z`;  chomp $tz;


#  Set up for the machine we're on.


if      (-d "/gryphon") {
    $site    = "gryphon";
    $wrkdir  = "/data/projects/phillippy/scratch/NIGHTLY";
    $gitrepo = "/data/projects/phillippy/scratch/NIGHTLY/canu";
}

elsif (-d "/data/walenzbp/NIGHTLY/") {
    $site    = "BPWI";
    $wrkdir  = "/data/walenzbp/NIGHTLY";
    $gitrepo = "/data/walenzbp/NIGHTLY/canu";
}

elsif (-d "/assembly/NIGHTLY/") {
    $site    = "BPWI";
    $wrkdir  = "/assembly/NIGHTLY";
    $gitrepo = "/assembly/NIGHTLY/canu";
}

else {
    die "Unknown site configuration.\n";
}


#  Now do something.


if ($oper eq "fetch") {
    fetch();
    exit(0);
}


if ($oper eq "checkout") {
    checkout(@ARGV);
    exit(0);
}


if ($oper eq "build") {
    build(@ARGV);
    exit(0);
}


if ($oper eq "assemble") {
    #  if canu.merge.out is "Already up-to-date." we can skip running stuff.
    assemble(@ARGV);
    exit(0);
}


if ($oper eq "submit") {
    submit(@ARGV);
    exit(0);
}


die "$0: unknown action '$oper'\n";
exit(0);






sub parseDate ($) {
    my $pathdate = shift @_;
    my $thisdate;
    my $lastdate;

    #  Until our Makefile is fixed, the ":" cannot be used, and so we need
    #  to adopt a slightly goofy timestamp name.

    if (defined($pathdate)) {
        if ($pathdate =~ m/^\d\d\d\d-\d\d-\d\d-\d\d\d\d$/) {
            $thisdate = $pathdate;
        } else {
            die "Malformed ddir '$pathdate'\n";
        }
    } else {
        my @v = localtime;
        $v[5] += 1900;
        $v[4]++;

        $v[5] = substr("0000$v[5]", -4);
        $v[4] = substr("0000$v[4]", -2);
        $v[3] = substr("0000$v[3]", -2);
        $v[2] = substr("0000$v[2]", -2);
        $v[1] = substr("0000$v[1]", -2);

        $thisdate = "$v[5]-$v[4]-$v[3]-$v[2]$v[1]";
    }

    open(F, "ls $wrkdir |");
    while (<F>) {
        chomp;
        if (m/^\d\d\d\d-\d\d-\d\d-\d\d\d\d$/) {
            if ($_ lt $thisdate) {
                $lastdate = $_;
            }
        }
    }
    close(F);

    #print STDERR "Working on '$thisdate' -- last found is '$lastdate'.\n";

    return($thisdate, $lastdate);
}



sub gitDate ($) {
    my $date = shift @_;

    if ((!defined($date)) || ($date eq "")) {
        return($date);
    }

    if ($date =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
        $date = "$1-$2-$3T$4:$5$tz";
    } else {
        die "Malformed date '$date' to gitDate().\n";
    }

    return($date);
}



sub fetch () {

    print STDERR "\n";
    print STDERR "SANITY FETCH\n";
    print STDERR "\n";

    if (! -d "$gitrepo") {
        print STDERR "  Cloning new repo '$gitrepo'.\n";
        system("mkdir -p $gitrepo")  if (! -d "$gitrepo");
        system("cd $gitrepo/.. ; git clone git\@github.com:marbl/canu.git > canu.clone.out 2>&1");
    }
    else {
        print STDERR "  Updating existing repo '$gitrepo'.\n";
        system("cd $gitrepo    ; git fetch > ../canu.fetch.out 2>&1");
        system("cd $gitrepo    ; git merge > ../canu.merge.out 2>&1");
    }
}



sub checkout (@) {
    my ($thisdate, $lastdate) = parseDate($_[0]);

    my $ldate = gitDate($lastdate);
    my $tdate = gitDate($thisdate);

    print STDERR "\n";
    print STDERR "SANITY CHECKOUT $wrkdir/$thisdate\n";
    print STDERR "\n";

    if (-d "$wrkdir/$thisdate/canu") {
        print STDERR "  $wrkdir/$thisdate/canu already exists.  Please remove to rerun.\n";
        return;
    }

    system("mkdir -p $wrkdir/$thisdate/canu") if (! -d "$wrkdir/$thisdate/canu");

    #  Clone the clone.

    system("cd $wrkdir/$thisdate/canu && rsync -a $gitrepo/ .");

    #  Find a git hash to grab the bits we want.

    my $hash = `cd $wrkdir/$thisdate/canu && git rev-list -n 1 --first-parent --before=\"$tdate\" master`;  chomp $hash;

    #  Checkout that hash.

    system("cd $wrkdir/$thisdate/canu && git checkout -b SANITY-$thisdate $hash");

    #  git clone $gitrepo $wrkdir/$thisdate/canu

    if ($ldate ne "") {
        print sTDERR "log\n";
        system("cd $wrkdir/$thisdate/canu && git log --after=\"$ldate\" --until=\"$tdate\" > $wrkdir/$thisdate/canu.updates");
    }

    print STDERR "  $thisdate checked out!\n";
}



sub build (@) {
    my ($thisdate, $lastdate) = parseDate($_[0]);

    print STDERR "\n";
    print STDERR "SANITY BUILD $wrkdir/$thisdate\n";
    print STDERR "\n";

    if (-e "$wrkdir/$thisdate/canu/src/make.err") {
        print STDERR "  $wrkdir/$thisdate/canu was already built once (successuflly or not).  Please cleanup first.\n";
        return;
    }

    system("cd $wrkdir/$thisdate/canu/src && gmake -j 12 > make.out 2> make.err");
}



sub assemble (@) {
    my ($thisdate, $lastdate) = parseDate(shift @_);

    print STDERR "\n";
    print STDERR "ASSEMBLE\n";
    print STDERR "\n";

    my $syst = `uname -s`;    chomp $syst;  #  OS implementation
    my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
    my $name = `uname -n`;    chomp $name;  #  Name of the system

    $arch = "amd64"  if ($arch eq "x86_64");
    $arch = "ppc"    if ($arch eq "Power Macintosh");

    my @spec;
    my $addresses;

    foreach my $arg (@_) {
        if ($arg =~ m/\@/) {
            $addresses  =   $arg    if (!defined($addresses));
            $addresses .= ",$arg"   if ( defined($addresses));
        }
        else {
            $arg = "$ENV{PWD}/$arg" if ($arg !~ m/^\//);
            push @spec, $arg;
        }
    }

    #open(F, "> $wrkdir/$thisdate/asm-done.sh");
    #print F "#!/bin/sh\n";
    #print F "\n";
    #print F "perl $wrkdir/sanity-asm-done.pl $wrkdir \$1 $thisdate\n";
    #close(F);

    #open(F, "> $wrkdir/$thisdate/all-done.sh");
    #print F "#!/bin/sh\n";
    #print F "\n";
    #print F "perl $wrkdir/sanity-all-done.pl $wrkdir $thisdate \$1 $addresses \$2 \\\n";
    #print F "| \\\n";
    #print F "tee \"$wrkdir/$thisdate/sanity-all-done.\$1\" \\\n";
    #print F "| \\\n";
    #print F "/usr/sbin/sendmail -i -t -f thebri\@gmail.com\n"  if ($site eq "JCVI");
    #print F "/usr/local/sbin/ssmtp thebri\@gmail.com\n"        if ($site eq "BPWI");
    #close(F);

    foreach my $s (@spec) {
        my @c = split '/', $s;
        my $n = $c[scalar(@c) - 1];

        $n =~ s/.spec$//;
        $n =~ s/.specFile$//;

        print STDERR "Submitting assembly '$n' for spec '$s'.\n";

        print STDERR "cd $wrkdir/$thisdate \\\n";
        print STDERR "&& \\\n";
        print STDERR "$wrkdir/$thisdate/canu/$syst-$arch/bin/canu -p $n -d $n -s $s\n";

        system("cd $wrkdir/$thisdate && $wrkdir/$thisdate/canu/$syst-$arch/bin/canu -p $n -d $n -s $s");
    }
}



sub submit (@) {
    my ($thisdate, $lastdate) = parseDate(shift @ARGV);
    my $seconds;

    if ($thisdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
        $seconds = timelocal(0, $5, $4, $3, $2 - 1, $1);
    }

    $seconds += 604800;  # one week
    $seconds += 86400;   # one day
    $seconds += 14400;   # four hours
    $seconds += 7200;    # two hours
    $seconds += 21600;   # six hours
    $seconds += 43200;   # twelve hours

    my @v = localtime($seconds);

    $v[5] += 1900;
    $v[4]++;

    $v[5] = substr("0000$v[5]", -4);
    $v[4] = substr("0000$v[4]", -2);
    $v[3] = substr("0000$v[3]", -2);
    $v[2] = substr("0000$v[2]", -2);
    $v[1] = substr("0000$v[1]", -2);

    #$v[2] = "00";
    #$v[1] = "01";

    my $nextdate = "$v[5]-$v[4]-$v[3]-$v[2]$v[1]";
    my $nexthold = "$v[5]$v[4]$v[3]$v[2]$v[1].00";
    my $nextname = "$v[4]$v[3]-$v[2]$v[1]";

    print STDERR "Submit next at date='$nextdate' hold='$nexthold' name='$nextname'\n";

    #system("qsub -cwd -j y -o $nextdate.err -A assembly-nightly -N snty$nextname -a $nexthold -b n sanity.sh $nextdate grid");
}


