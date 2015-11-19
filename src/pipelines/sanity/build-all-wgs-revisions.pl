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

#--------------------------------------------------------------------------------

use strict;

my $wgssvn  = "/work/NIGHTLY/SVN-wgs";
my $kmersvn = "/work/NIGHTLY/SVN-kmer";

my $wgsbyd  = "/work/NIGHTLY/bydate";

system("mkdir $wgsbyd")  if (! -d "$wgsbyd");

my %revToDate;

if (-d "$wgsbyd/TIP") {
    print STDERR "Updating to latest.  DON'T FORGET TO sanity.pl rsync FIRST!\n";
    system("cd $wgsbyd/TIP/src ; svn update > update.err");
} else {
    print STDERR "Checking out latest.\n";
    system("mkdir $wgsbyd/TIP");
    system("cd $wgsbyd/TIP ; svn co file://$wgssvn/trunk/src src > checkout.err");
}

print STDERR "Reading CA commit logs.\n";

open(F, "cd $wgsbyd/TIP/src ; svn log |");
while (<F>) {
    #print STDERR $_;

    #  CVS relic
    #if (m/^date:\s+(....)\/(..)\/(..)\s+(..):(..):(..);\s+author:\s+(.*);\s+state/) {
    #    my $fildate = "$1-$2-$3-$4$5";
    #    my $cvsdate = "$1/$2/$3 $4:$5:$6 GMT";
    #    my $svndate = "{$1$2$3T$4$5Z}";
    #    $logdate{"$fildate\0$cvsdate\0$svndate"}++;
    #}

    if (m/^(r\d+)\s+\|\s+(.+)\s+\|\s+(\d\d\d\d)-(\d\d)-(\d\d)\s(\d\d):(\d\d):(\d\d)\s(.\d\d\d\d)\s/) {
        my $fildate = "$3-$4-$5-$6$7";
        my $svndate = "{$3$4$5T$6$7Z}";

        $revToDate{$1} = "$fildate\0$svndate\0$2";
    }
}
close(F);



my @keys = sort keys %revToDate;

print STDERR "Checking out versions.\n";

foreach my $wgsrev (@keys) {
    my ($dirdate, $svndate, $author) = split '\0', $revToDate{$wgsrev};

    #next if (-e "$wgsbyd/$dirdate");
    #next if (-e "$wgsbyd/$dirdate.tar");
    #next if (-e "$wgsbyd/$dirdate.tar.gz");
    #next if (-e "$wgsbyd/$dirdate.tar.bz2");
    #next if (-e "$wgsbyd/$dirdate.tar.xz");

    next if ($dirdate eq "2004-04-14-1349");  #  These are all one
    next if ($dirdate eq "2004-04-14-1350");  #  slow commit, which
    next if ($dirdate eq "2004-04-14-1351");  #  we capture in
    next if ($dirdate eq "2004-04-14-1352");  #  2004-04-14-1354
    next if ($dirdate eq "2004-04-14-1353");  #

    next if ($dirdate le "2010-00");

    next if ($dirdate eq "2011-09-06-1707");  #  Doesn't compile

    next if ($dirdate eq "2012-02-01-2010");  #  Broken GWS overlapStore
    next if ($dirdate eq "2012-02-01-2012");  #  Broken GWS overlapStore
    next if ($dirdate eq "2012-02-01-2015");  #  Broken GWS overlapStore

    next if ($dirdate le "2014-00");

    system("mkdir -p $wgsbyd/$dirdate")       if (! -e "$wgsbyd/$dirdate");
    system("ln -s $dirdate $wgsbyd/$wgsrev")  if (! -e "$wgsbyd/$wgsrev");

    #  Checkout the assembler

    if (! -e "$wgsbyd/$dirdate/src") {
        print "$dirdate -- $wgsrev CHECKOUT\n";

        system("cd $wgsbyd/$dirdate ; svn co -r $wgsrev file://$wgssvn/trunk/src src > src.$wgsrev.err");
    } else {
        print "$dirdate -- $wgsrev UPDATE\n";

        system("cd $wgsbyd/$dirdate/src ; svn update -r $wgsrev");
    }

    #  Checkout kmer, or, more probably, just link to one that exists already

    my $kmerev = "0000";

    open(L, "svn info -r \"$svndate\" file://$kmersvn/trunk |") or die "Failed to get info: $!\n";
    while (<L>) {
        $kmerev = $1  if (m/^Revision:\s+(\d+)/);
    }
    close(L);

    die "Didn't find a kmer revision in $wgsbyd/$dirdate/kmer/kmer.err.\n"  if ($kmerev eq "0000");

    if (! -e "$wgsbyd/kmer$kmerev") {
        print "$dirdate -- kmer$kmerev CHECKOUT\n";

        system("cd $wgsbyd/$dirdate && svn co -r \"$svndate\" file://$kmersvn/trunk $wgsbyd/kmer$kmerev > $wgsbyd/kmer$kmerev.err 2>&1");
    }

    system("ln -s ../kmer$kmerev $wgsbyd/$dirdate/kmer")   if (! -e "$wgsbyd/$dirdate/kmer");

    #  Compile

    if (! -e "$wgsbyd/kmer$kmerev/build.err") {
        print "$dirdate -- kmer$kmerev COMPILE\n";

        if (! -e "$wgsbyd/kmer$kmerev/Makefile") {
            system("cd $wgsbyd/kmer$kmerev && sh configure.sh > configure.err 2>&1");
        }

        open(F, "> $wgsbyd/kmer$kmerev/build.err");
        print F "Waiting to compile.\n";
        close(F);

        system("qsub -N kmer$kmerev -cwd -j y -o $wgsbyd/kmer$kmerev/build.err -q vomit.q -wd $wgsbyd/kmer$kmerev -b y gmake install");
        #system("cd $wgsbyd/kmer$kmerev && gmake install > build.err 2>&1 &");
    }

    if (! -e "$wgsbyd/$dirdate/src/build.err") {
        print "$dirdate -- $wgsrev COMPILE\n";

        open(F, "> $wgsbyd/$dirdate/src/build.err");
        print F "Waiting to compile.\n";
        close(F);

        system("qsub -N wgs$wgsrev -hold_jid kmer$kmerev -cwd -j y -o $wgsbyd/$dirdate/src/build.err -q vomit.q -wd $wgsbyd/$dirdate/src -b y gmake");
        #system("cd $wgsbyd/$dirdate/src && gmake > build.err 2>&1");
    }

    #  Archive

    #if (! -e "$wgsbyd/$dirdate.tar.bz2") {
    #    Use sge, hold on both $kmerev and $wgsrev
    #    system("cd $wgsbyd && tar -cf - $dirdate | xz -9vc > $dirdate.tar.bz2 && mv $wgsbyd/$dirdate $wgsbyd/$dirdate.DELETE &");
    #}
}

