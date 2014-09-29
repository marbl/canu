#!/usr/bin/perl
#
#  Checkout and tar up src for each commit.
#
#  The first half comes from get-commit-logs.pl.
#
#--------------------------------------------------------------------------------

use strict;

my $wgscvs  = "/work/NIGHTLY/wgs-assembler-cvs";
my $kmersvn = "/work/NIGHTLY/kmer-svn";

my $wgsbyd  = "/work/NIGHTLY/bydate";

if (! -d "$wgsbyd") {
    system("mkdir $wgsbyd");
}

my %logdate;

if (-d "$wgsbyd/TIP") {
    print STDERR "Updating to latest.  DON'T FORGET TO sanity.pl rsync FIRST!\n";
    system("cd $wgsbyd/TIP/src ; cvs update > update.err");
} else {
    print STDERR "Checking out latest.\n";
    system("mkdir $wgsbyd/TIP");
    system("cd $wgsbyd/TIP ; cvs -d $wgscvs co src > src.err");
}

print STDERR "Reading commit logs.\n";

open(F, "cd $wgsbyd/TIP/src ; cvs log -N -S |");
while (<F>) {
    #print STDERR $_;
    if (m/^date:\s+(....)\/(..)\/(..)\s+(..):(..):(..);\s+author:\s+(.*);\s+state/) {
        my $fildate = "$1-$2-$3-$4$5";
        my $cvsdate = "$1/$2/$3 $4:$5:$6 GMT";
        my $svndate = "{$1$2$3T$4$5Z}";
        #print "$fildate -- $cvsdate\n";
        $logdate{"$fildate\0$cvsdate\0$svndate"}++;
    }
}
close(F);

my @keys = sort keys %logdate;

print STDERR "Checking out versions.\n";

foreach my $d (@keys) {
    my ($f, $c, $s) = split '\0', $d;

    next if (-e "$wgsbyd/$f.tar");
    next if (-e "$wgsbyd/$f.tar.gz");
    next if (-e "$wgsbyd/$f.tar.bz2");

    next if ($f eq "2004-04-14-1349");  #  These are all one
    next if ($f eq "2004-04-14-1350");  #  slow commit, which
    next if ($f eq "2004-04-14-1351");  #  we capture in
    next if ($f eq "2004-04-14-1352");  #  2004-04-14-1354
    next if ($f eq "2004-04-14-1353");  #

    next if ($f le "2010-00");

    next if ($f eq "2011-09-06-1707");  #  Doesn't compile

    next if ($f eq "2012-02-01-2010");  #  Broken GWS overlapStore
    next if ($f eq "2012-02-01-2012");  #  Broken GWS overlapStore
    next if ($f eq "2012-02-01-2015");  #  Broken GWS overlapStore

    if (! -e "$wgsbyd/$f") {
        system("mkdir -p $wgsbyd/$f");
    }

    #  Checkout the assembler

    if (! -e "$wgsbyd/$f/src") {
        print "$f -- $c -- checkout wgs\n";

        system("cd $wgsbyd/$f ; cvs -r -d $wgscvs  -z3 co -N -D \"$c\" src > src.err");
    } else {
        print "$f -- $c -- update wgs\n";

        system("cd $wgsbyd/$f ; cvs update");
    }

    #  Checkout kmer

    if (! -e "$wgsbyd/$f/kmer") {
        print "$f -- $c -- checkout kmer\n";

        system("cd $wgsbyd/$f && svn co  -r \"$s\" file://$kmersvn/trunk kmer > kmer.err 2>&1");

        #  Analyze the kmer checkout log to decide which revision we have.  Create a link to it.

        my $revision = "0000";

        open(L, "< $wgsbyd/$f/kmer.err") or die "Failed to open '$wgsbyd/$f/kmer.err'\n";
        while (<L>) {
            chomp;
            
            if (m/^Checked\s+out\s+revision\s+(\d+).$/) {
                $revision = $1;
            }
        }
        close(L);
 
        if ($revision eq "0000") {
            die "Didn't find a kmer revision in $wgsbyd/$f/kmer/kmer.err.\n";

        } elsif (! -e "$wgsbyd/kmer$revision") {
            print "$f -- $c -- kmer$revision is NEW\n";

            system("mv $wgsbyd/$f/kmer     $wgsbyd/kmer$revision");
            system("mv $wgsbyd/$f/kmer.err $wgsbyd/kmer$revision.err");

            system("ln -s $wgsbyd/kmer$revision $wgsbyd/$f/kmer");

        } else {
            print "$f -- $c -- kmer$revision\n";

            system("mv $wgsbyd/$f/kmer     $wgsbyd/kmer$revision.$f.DELETE");
            system("mv $wgsbyd/$f/kmer.err $wgsbyd/kmer$revision.$f.DELETE.err");

            system("ln -s $wgsbyd/kmer$revision $wgsbyd/$f/kmer");
        }

        if (! -e "$wgsbyd/kmer$revision/FreeBSD-amd64/bin/meryl") {
            print "$f -- $c -- kmer$revision -- compile\n";

            if (! -e "$wgsbyd/kmer$revision/Makefile") {
                system("cd $wgsbyd/kmer$revision && sh configure.sh > configure.err 2>&1");
            }

            system("cd $wgsbyd/kmer$revision && gmake install > build.err 2>&1");
        }
    }

    #  Compile


    if (! -e "$wgsbyd/$f/FreeBSD-amd64/bin/gatekeeper") {
        print "$f -- $c -- compile wgs\n";

        #system("cd $wgsbyd/$f/src && gmake > build.err 2>&1");
        system("qsub -cwd -j y -o $wgsbyd/$f/src/build.err -q vomit.q -wd $wgsbyd/$f/src -b y gmake");
    }


    #if (! -e "$wgsbyd/$f.tar.bz2") {
    #    system("cd $wgsbyd && tar -cf - $f | bzip2 -9vc > $f.tar.bz2 && mv $wgsbyd/$f $wgsbyd/$f.DELETE &");
    #}
}

