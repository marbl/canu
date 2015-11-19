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

my $bgn = shift @ARGV;
my $end = shift @ARGV;

my @dates;

open(F, "ls -r *bz2 |");
while (<F>) {
    chomp;

    next  if (defined($bgn) && ($_ lt $bgn));
    next  if (defined($end) && ($end lt $_));

    push @dates, $_;
}

$bgn = "START_OF_TIME"  if (!defined($bgn));
$end = "END_OF_TIME"    if (!defined($end));

print STDERR "Building ", scalar(@dates), " revisions from $bgn to $end.\n";

while (scalar(@dates) > 0) {
    my $bz2 = shift @dates;
    my $dir;

    if ($bz2 =~ m/(\d\d\d\d-\d\d-\d\d-\d\d\d\d).tar.bz2/) {
        $dir = $1;
    } else {
        die "Nope bz2='$bz2'\n";
    }

    if (! -e "$dir") {
        print STDERR "Uncompressing $bz2\n";
        system("bzip2 -dc $bz2 | tar -xf -") and die "Failed to uncompress '$bz2'\n";
        next;
    }

    #if (! -e "$dir/kmer/FreeBSD-amd64") {
    #    print STDERR "Building $dir/kmer\n";
    #    if ($dir le "2000-00-00-0000") {
    #        die;
    #    } else {
    #        system("cd $dir/kmer && gmake install > errs 2>&1") and die "Failed to build '$dir/kmer'\n";
    #    }
    #}

    next if ($dir le "2008-10-10-0304");  #  ALL the early assemblers fail to build; no kmer

    next if ($dir eq "2008-10-29-1053");  #  buildRefUnitigs is not committed
    next if ($dir eq "2008-10-29-1619");  #  buildRefUnitigs is not committed
    next if ($dir eq "2008-10-29-1649");  #  buildRefUnitigs is not committed
    next if ($dir eq "2008-10-29-1651");  #  buildRefUnitigs is not committed
    next if ($dir eq "2008-10-29-1719");  #  buildRefUnitigs is not committed
    next if ($dir eq "2008-10-29-1720");  #  buildRefUnitigs is not committed
    next if ($dir eq "2008-10-30-0449");  #  buildRefUnitigs is not committed
    next if ($dir eq "2008-11-02-0628");
    next if ($dir eq "2009-01-16-1639");
    next if ($dir eq "2009-01-16-1643");
    next if ($dir eq "2009-01-16-1646");
    next if ($dir eq "2009-01-16-1647");
    next if ($dir eq "2009-01-16-1657");
    next if ($dir eq "2009-03-06-2012");
    next if ($dir eq "2009-03-31-2004");
    next if ($dir eq "2009-03-31-2032");
    next if ($dir eq "2009-04-08-1725");
    next if ($dir eq "2009-04-09-1401");
    next if ($dir eq "2009-06-10-1805");
    next if ($dir eq "2009-06-24-1205");
    next if ($dir eq "2009-07-06-2003");
    next if ($dir eq "2009-07-27-0806");
    next if ($dir eq "2009-07-28-1230");
    next if ($dir eq "2009-08-04-1103");
    next if ($dir eq "2009-08-04-1105");
    next if ($dir eq "2009-08-06-1136");
    next if ($dir eq "2009-08-14-1103");  #  BREAKS WITH rename-to-c++
    next if ($dir eq "2009-08-28-17597");
    next if ($dir eq "2009-09-04-2024");
    next if ($dir eq "2009-09-04-2025");
    next if ($dir eq "2009-09-07-0740");
    next if ($dir eq "2009-09-09-0340");
    next if ($dir eq "2009-09-09-0745");
    next if ($dir eq "2009-10-12-0619");
    next if ($dir eq "2009-10-27-1244");
    next if ($dir eq "2009-11-19-1501");
    next if ($dir eq "2009-12-03-0119");
    next if ($dir eq "2009-12-19-0536");
    next if ($dir eq "2010-01-14-0044");
    next if ($dir eq "2010-01-26-0351");  #  AS_BOG_BestOverlapGraph() prototype mismatch
    next if ($dir eq "2010-02-04-2156");  #  BOG fails to build
    next if ($dir eq "2010-02-09-1812");  #  markUniqueUnique() bad header file
    next if ($dir eq "2010-02-17-0132");  #  pointer/object confusion
    next if ($dir eq "2010-03-22-2008");  #  AS_UTL_writeFastQ() not defined
    next if ($dir eq "2010-09-28-1055");  #  Makefile refers to missing overlapStore_genomeLength.c
    next if ($dir eq "2010-09-28-1720");  #  Makefile refers to missing overlapStore_genomeLength.c
    next if ($dir eq "2010-10-01-1343");
    next if ($dir eq "2010-10-04-0852");
    next if ($dir eq "2010-10-15-0246");  #  const char * vs char * in AS_BOG
    next if ($dir eq "2010-12-08-1242");  #  Makefile refers to deleted AS_CGW/AS_CGW_dump.c
    next if ($dir eq "2010-12-08-1243");  #  Makefile refers to deleted AS_CGW/AS_CGW_dump.c
    next if ($dir eq "2010-12-08-1245");  #  Makefile refers to deleted AS_CGW/AS_CGW_dump.c
    next if ($dir eq "2010-12-09-0405");  #  Makefile refers to deleted AS_CGW/AS_CGW_dump.c
    next if ($dir eq "2010-12-10-1328");  #  Makefile refers to deleted AS_CGW/AS_CGW_dump.c
    next if ($dir eq "2011-01-25-1356");
    next if ($dir eq "2011-02-11-0548");
    next if ($dir eq "2011-03-08-2118");  #  Undefined constant in AS_BOG
    next if ($dir eq "2011-06-17-1303");
    next if ($dir eq "2011-08-01-1835");
    next if ($dir eq "2011-08-01-2033");
    next if ($dir eq "2011-08-01-2100");
    next if ($dir eq "2011-08-01-2159");
    next if ($dir eq "2011-08-01-2233");
    next if ($dir eq "2011-08-02-0221");
    next if ($dir eq "2011-08-02-0223");
    next if ($dir eq "2011-08-02-0225");
    next if ($dir eq "2011-08-02-0318");
    next if ($dir eq "2011-08-04-1708");
    next if ($dir eq "2011-08-24-2114");
    next if ($dir eq "2011-08-25-0240");
    next if ($dir eq "2011-08-30-1332");
    next if ($dir eq "2011-08-30-1812");
    next if ($dir eq "2011-08-30-1813");
    next if ($dir eq "2011-08-30-2309");
    next if ($dir eq "2011-09-02-1714");  #  ouch, broken through 2011-09-06-1707
    next if ($dir eq "2011-09-02-1825");  #  ouch.
    next if ($dir eq "2011-09-02-2204");  #  ouch.
    next if ($dir eq "2011-09-03-0129");  #  ouch.
    next if ($dir eq "2011-09-03-0408");  #  ouch.
    next if ($dir eq "2011-09-03-0736");  #  ouch.
    next if ($dir eq "2011-09-03-0813");  #  ouch.
    next if ($dir eq "2011-09-04-0101");  #  ouch.
    next if ($dir eq "2011-09-05-1649");  #  ouch.
    next if ($dir eq "2011-09-05-2123");  #  ouch.
    next if ($dir eq "2011-09-06-0111");  #  ouch.
    next if ($dir eq "2011-09-06-0215");  #  ouch.
    next if ($dir eq "2011-09-06-0947");  #  ouch.
    next if ($dir eq "2011-09-06-1415");  #  ouch.
    next if ($dir eq "2011-09-06-1422");  #  ouch.
    next if ($dir eq "2011-09-06-1506");  #  ouch.
    next if ($dir eq "2011-09-06-1524");  #  ouch.
    next if ($dir eq "2011-09-06-1527");  #  ouch.
    next if ($dir eq "2011-09-06-1601");  #  ouch.
    next if ($dir eq "2011-09-06-1641");  #  ouch.
    next if ($dir eq "2011-09-06-1707");  #  ouch.
    next if ($dir eq "2011-12-09-2302");  #  Assigning 64-bit allocation to 32-bit pointer in AS_BAT


    if (! -e "$dir/kmer/FreeBSD-amd64") {
        die "Didn't find kmer in $dir\n";
    }

    if (! -e "$dir/FreeBSD-amd64/bin/gatekeeper") {
        print STDERR "Building $dir/src\n";

        if      ($dir le "2009-02-03-2104") {
            #  These old builds REQUIRE that SITE_NAME be set explicitly.
            $ENV{'SITE_NAME'} = "LOCAL";
            system("cd $dir/src && gmake > errs 2>&1") and die "Failed to build '$dir/src'\n";
            undef $ENV{'SITE_NAME'};

        } elsif ($dir lt "2009-06-12-1741") {
            #  These builds are before the C++ hack
            system("cd $dir/src && gmake > errs 2>&1") and die "Failed to build '$dir/src'\n";

        } elsif ($dir le "2009-08-05-2205") {
            #  These builds need to rename C to C++ to build.
            #  The rename script doesn't go away until (after) 2009-09-30-1832
            system("cd $dir/src && sh rename-to-c++.sh") and die "Failed to rename '$dir/src'\n";
            system("cd $dir/src && gmake > errs 2>&1") and die "Failed to build '$dir/src'\n";

        } else {
            system("cd $dir/src && gmake > errs 2>&1") and die "Failed to build '$dir/src'\n";
        }
    }
}
