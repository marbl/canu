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

#  The only four globals configurables:
#
my $site    = undef;
my $wrkdir  = undef;
my $wgssvn  = undef;
my $kmersvn = undef;

if      (-d "/usr/local/projects/BIOINFO/ASSEMBLY/NIGHTLY") {
    $site    = "JCVI";
    $wrkdir  = "/usr/local/projects/BIOINFO/ASSEMBLY/NIGHTLY";
    $wgssvn  = "/usr/local/projects/BIOINFO/ASSEMBLY/NIGHTLY/SVN-wgs";
    $kmersvn = "/usr/local/projects/BIOINFO/ASSEMBLY/NIGHTLY/SVN-kmer";
} elsif (-d "/work/NIGHTLY/") {
    $site    = "BPWI";
    $wrkdir  = "/work/NIGHTLY";
    $wgssvn  = "/work/NIGHTLY/SVN-wgs";
    $kmersvn = "/work/NIGHTLY/SVN-kmer";
} else {
    die "Unknown site configuration.\n";
}


#  Command line options are
#
#    'oper'   -- thing to do.
#
#    'ddir' -- properly formatted directory name to do it in.  op ==
#    checkout creates this directory.  If still supplied, it will
#    force a checkout from that date.
#
#  rsync
#     update the local repositories in 'wgs-assembler-cvs' and 'kmer-svn'.
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
#  assemble date label specFile .... email ...
#     Launches runCA on each specFile, running it in directory 'date'.
#     Each specFile must have a unique name; the assembly is named
#     after the specFile.  The label is used only for reporting via email.
#     At the end, a diff is generated to POINTERS/$prefix.last and
#     POINTERS/$prefix.reference.  If the assembly finished, POINTERS/$prefix.last
#     is updated.
#

{
    my $oper   = shift @ARGV;
    my $ddir   = shift @ARGV;

    my ($thisdate, $lastdate) = parseDate($ddir);

    if ($oper eq "rsync") {
        system("mkdir -p $wgssvn")  if (! -d "$wgssvn");
        system("mkdir -p $kmersvn") if (! -d "$kmersvn");

        #system("cd $wgscvs  && rsync -av rsync://wgs-assembler.cvs.sourceforge.net/cvsroot/wgs-assembler/\* . > rsync.out 2>&1");
        system("cd $wgssvn  && rsync -av svn.code.sf.net::p/wgs-assembler/svn/ . > rsync.out 2>&1");
        system("cd $kmersvn && rsync -av svn.code.sf.net::p/kmer/code/         . > rsync.out 2>&1");

    } elsif ($oper eq "checkout") {
        checkoutAndLog($thisdate, $lastdate, "file://$kmersvn/trunk",    "kmer", "kmer");
        checkoutAndLog($thisdate, $lastdate, "file://$wgssvn/trunk/src", "src", "src");

    } elsif ($oper eq "build") {
        buildKmer($thisdate);
        buildCA($thisdate);

    } elsif ($oper eq "assemble") {
        assemble($thisdate, @ARGV);

    } else {
        die "$0: unknown action '$oper'\n";
    }

    exit(0);
}


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

    print STDERR "Working on '$thisdate' -- last found is '$lastdate'.\n";

    return($thisdate, $lastdate);
}


sub checkoutAndLog ($$) {
    my $thisdate = shift @_;
    my $lastdate = shift @_;
    my $repo     = shift @_;  #  Path to local repository
    my $target   = shift @_;  #  Thing to checkout
    my $path     = shift @_;  #  Where to put it (in wgs/)

    print STDERR "Checking out $wrkdir/$thisdate (repo=$repo path=$path)\n";

    if (-d "$wrkdir/$thisdate/wgs/$path") {
        print STDERR "$wrkdir/$thisdate/wgs/$path already exists.  Please remove to rerun.\n";
        return;
    }

    system("mkdir -p $wrkdir/$thisdate/wgs") if (! -d "$wrkdir/$thisdate/wgs");

    #  Convert time-names to dates that svn can use.

    my $thisdatesvn;
    my $lastdatesvn;

    #  NOT tested with SVN.  The old format was "-04:00"; new format is "-0400".
    my $tz = `date +%z`;  chomp $tz;

    if ($thisdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
        $thisdatesvn = "$1-$2-$3T$4:$5$tz";
    } else {
    }

    if (defined($lastdate)) {
        if ($lastdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
            $lastdatesvn = "$1-$2-$3T$4:$5$tz";
        } else {
        }
    }

    system("cd $wrkdir/$thisdate/wgs && svn co  -r \"{$thisdatesvn}\" $repo $target > $path.checkout.err 2>&1");

    #  This is annoying.  SVN log will report changes inclusive to revisions.  -r 5:9 will report
    #  changes made in revisions 5 through 9.  When you give it a date, it finds the revision that
    #  was active on that date.  What we want here, though, is the changes SINCE that date (or,
    #  since revision 5, up to revision 9).
    #
    #  To get around this, we first get the logs, but scan for the lowest and highest revision
    #  numbers, then dump correctly.
    #
    if ($lastdate ne "") {
        my $loRev;
        my $hiRev;

        print "svn log -v $repo -r \"{$lastdatesvn}:{$thisdatesvn}\"\n";
        open(F, "cd $wrkdir/$thisdate/wgs && svn log -v $repo -r \"{$lastdatesvn}:{$thisdatesvn}\" |");
        while (<F>) {
            if (m/^r(\d+)\s+\|\s+/) {
                $loRev = $1  if ((!defined($loRev)) || ($1 < $loRev));
                $hiRev = $1  if ((!defined($hiRev)) || ($hiRev < $1));
            }
        }
        close(F);

        $loRev++  if (defined($loRev));

        print STDERR "loRev='$loRev' hiRev='$hiRev'\n";

        if (defined($loRev) && defined($hiRev) && ($loRev < $hiRev)) {
            print "svn log -v $repo -r $loRev:$hiRev\n";
            system("cd $wrkdir/$thisdate/wgs && svn log -v $repo -r $loRev:$hiRev > $path.updates");
        }
    }

    print STDERR "$thisdate checked out!\n";
}



sub buildKmer ($) {
    my $thisdate = shift @_;

    print STDERR "Building KMER $wrkdir/$thisdate\n";

    if (-e "$wrkdir/$thisdate/wgs/kmer/make.err") {
        print STDERR "$wrkdir/$thisdate/wgs/kmer was already built once (successuflly or not).  Please cleanup first.\n";
    } else {
        system("cd $wrkdir/$thisdate/wgs/kmer && sh configure.sh > configure.out 2> configure.err");
        system("cd $wrkdir/$thisdate/wgs/kmer && gmake install   > make.out.raw  2> make.err.raw");

        my %lines;

        open(F, "< $wrkdir/$thisdate/wgs/kmer/make.err.raw");
        open(G, "> $wrkdir/$thisdate/wgs/kmer/make.err");
        while (<F>) {
            chomp;
            next if (m/^ar:\s+creating/);
            if (!exists($lines{$_})) {
                $lines{$_}++;
                print G "$_\n";
            }
        }
        close(F);
        close(G);
    }
}


sub buildCA ($) {
    my $thisdate = shift @_;

    print STDERR "Building CA $wrkdir/$thisdate\n";

    if (-e "$wrkdir/$thisdate/wgs/src/make.err") {
        print STDERR "$wrkdir/$thisdate/wgs/src was already built once (successuflly or not).  Please cleanup first.\n";
    } else {
        #  Temporary hack to handle the C to C++ renaming.
        #  This was broken on "2009/06/10 17:34:22".
        #  This was fixed  on "2009/08/06 11:36:54"
        #
        if ((-e "$wrkdir/$thisdate/wgs/src/rename-to-c++.sh") &&
            (-s "$wrkdir/$thisdate/wgs/src/c_make.gen" == 14877)) {
            system("cd $wrkdir/$thisdate/wgs/src && sh rename-to-c++.sh");
        }

        system("cd $wrkdir/$thisdate/wgs/src && gmake > make.out.raw 2> make.err.raw");

        my %lines;

        open(F, "< $wrkdir/$thisdate/wgs/src/make.err.raw");
        open(G, "> $wrkdir/$thisdate/wgs/src/make.err");
        while (<F>) {
            chomp;
            next if (m/^ar:\s+creating/);
            if (!exists($lines{$_})) {
                $lines{$_}++;
                print G "$_\n";
            }
        }
        close(F);
        close(G);
    }
}



sub assemble ($$@) {
    my $thisdate  = shift @_;
    my $label     = shift @_;
    my $holds_asms;
    my $names_asms;
    my $holds_done;

    my $arch;

    #  Ripped from runCA/util.pl
    my $syst = `uname -s`;    chomp $syst;  #  OS implementation
    my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
    my $name = `uname -n`;    chomp $name;  #  Name of the system

    $arch = "amd64"  if ($arch eq "x86_64");
    $arch = "ppc"    if ($arch eq "Power Macintosh");
    #  end of rip.


    #  Figure out which are specfiles and which are email addresses.
    #
    my @spec;
    my @addresses;
    my $addresses;

    foreach my $arg (@_) {
        if ($arg =~ m/\@/) {
            push @addresses, $arg;
        } else {
            $arg = "$ENV{PWD}/$arg" if ($arg !~ m/^\//);
            push @spec, $arg;
        }
    }

    $addresses = join ',', @addresses;

    open(F, "> $wrkdir/$thisdate/asm-done.sh");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "perl $wrkdir/sanity-asm-done.pl $wrkdir \$1 $thisdate\n";
    close(F);

    open(F, "> $wrkdir/$thisdate/all-done.sh");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "perl $wrkdir/sanity-all-done.pl $wrkdir $thisdate \$1 $addresses \$2 \\\n";
    print F "| \\\n";
    print F "tee \"$wrkdir/$thisdate/sanity-all-done.\$1\" \\\n";
    print F "| \\\n";
    print F "/usr/sbin/sendmail -i -t -f celera_assembler_test\@jcvi.org\n" if ($site eq "JCVI");
    print F "/usr/local/sbin/ssmtp thebri\@gmail.com\n" if ($site eq "BPWI");
    close(F);

    foreach my $s (@spec) {
        my $n;

        {
            my @c = split '/', $s;
            $n = $c[scalar(@c) - 1];
            $n =~ s/.spec$//;
            $n =~ s/.specFile$//;
        }


        print STDERR "----------------------------------------\n";
        print STDERR "Submitting assembly '$n'.\n";

        system("mkdir $wrkdir/$thisdate/$n") if (! -d "$wrkdir/$thisdate/$n");

        my $jl = "CAini_${n}_$$";  #  Name of the launcher
        my $jn = "CAtst_${n}_$$";  #  Name of the asm-done

        open(F, "> $wrkdir/$thisdate/$n/launch-assembly.sh");
        print F "#!/bin/sh\n";
        print F "\n";
        print F "#  runCA checks if this is set to decide if it is on the grid or not.  We want\n";
        print F "#  to pretend we are NOT on the grid, so runCA will resubmit itself immediately.\n";
        print F "#\n";
        print F "unset SGE_TASK_ID\n";
        print F "\n";
        print F "\n";
        print F "#  Attempt to (re)configure SGE.  For reasons Bri doesn't know,\n";
        print F "#  jobs submitted to SGE, and running under SGE, fail to read his\n";
        print F "#  .tcshrc (or .bashrc, limited testing), and so they don't setup\n";
        print F "#  SGE (or ANY other paths, etc) properly.  For the record,\n";
        print F "#  interactive SGE logins (qlogin, etc) DO set the environment.\n";
        print F "\n";
        print F ". \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n";
        print F "\n";
        print F "\n";
        print F "#  Submit runCA to the grid.  Do not set any runCA parameters here; they\n";
        print F "#  will override ALL specFile parameters.  scriptOnGrid must be set, otherwise\n";
        print F "#  the assembly will be performed with this call, instead of just launched\n";
        print F "#  to the grid.\n";
        print F "\n";
        print F "perl $wrkdir/$thisdate/wgs/$syst-$arch/bin/runCA \\\n";
        print F "  sgePropagateHold=$jn scriptOnGrid=1 \\\n";
        print F "  -p $n -d $wrkdir/$thisdate/$n -s $s\n";
        print F "\n";
        print F "#  Once that runCA finishes, we've updated the hold on $jn, and so can release\n";
        print F "#  our user hold on it.\n";
        print F "#\n";
        print F "qrls -h u $jn\n";
        print F "\n";
        print F "\n";
        close(F);

        #  A separate script is used to launch runCA, so that we can get a user hold on it.  This is messy.
        #  1)  Submit the launcher, holding it.
        #  2)  Submit the asm-done, holding it too.  Make it also hold_jid on the launcher.
        #  3)  Release the launcher.
        #  4)  The launcher submits runCA, with sgePropagateHold.  When that runCA finishes, the assembly
        #      is now running, AND it has updated the hold_jid on asm-done.
        #  5)  We can now release the asm-done.
        #
        #  Steps 1 and 2 are done here, step 3 is at the very end, steps 4 and 5 are done in the launch-assembly.sh above.

        if      ($site eq "JCVI") {
            system("cd $wrkdir/$thisdate/$n && qsub -P 334007 -A CAsanity -b n -cwd -j y -o $wrkdir/$thisdate/$n/launch-assembly.err -l fast -h               -N $jl $wrkdir/$thisdate/$n/launch-assembly.sh");
            system("cd $wrkdir/$thisdate/$n && qsub -P 334007 -A CAsanity -b n -cwd -j y -o $wrkdir/$thisdate/$n/asm-done.err        -l fast -h -hold_jid $jl -N $jn $wrkdir/$thisdate/asm-done.sh $n");
        } elsif ($site eq "BPWI") {
            system("cd $wrkdir/$thisdate/$n && qsub           -A CAsanity -b n -cwd -j y -o $wrkdir/$thisdate/$n/launch-assembly.err         -h               -N $jl $wrkdir/$thisdate/$n/launch-assembly.sh");
            system("cd $wrkdir/$thisdate/$n && qsub           -A CAsanity -b n -cwd -j y -o $wrkdir/$thisdate/$n/asm-done.err                -h -hold_jid $jl -N $jn $wrkdir/$thisdate/asm-done.sh $n");
        }

        if (defined($holds_done)) {
            $holds_done .= ",$jn";
        } else {
            $holds_done  = "$jn";
        }

        if (defined($holds_asms)) {
            $holds_asms .= ",$jl";
            $names_asms .= ",$n";
        } else {
            $holds_asms  = "$jl";
            $names_asms .= "$n";
        }
    }

    if      ($site eq "JCVI") {
        system("cd $wrkdir/$thisdate && qsub -P 334007 -A CAsanity -b n -cwd -j y -o $wrkdir/$thisdate/all-done.err -l fast -hold_jid $holds_done -N CAfin_$thisdate $wrkdir/$thisdate/all-done.sh $label $names_asms");
    } elsif ($site eq "BPWI") {
        system("cd $wrkdir/$thisdate && qsub           -A CAsanity -b n -cwd -j y -o $wrkdir/$thisdate/all-done.err         -hold_jid $holds_done -N CAfin_$thisdate $wrkdir/$thisdate/all-done.sh $label $names_asms");
    }

    #  Now, release everything to run.

    system("qrls -h u $holds_asms");
}

exit(0);
