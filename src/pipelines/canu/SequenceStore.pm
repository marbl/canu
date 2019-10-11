
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
 #    Brian P. Walenz beginning on 2018-APR-23
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::SequenceStore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getNumberOfReadsInStore
             getNumberOfBasesInStore
             getSizeOfSequenceStore
             getExpectedCoverage
             generateReadLengthHistogram
             checkSequenceStore);

use strict;
use warnings "all";
no  warnings "uninitialized";

use Cwd qw(getcwd);

use canu::Defaults;
use canu::Execution;

use canu::Report;
use canu::Output;

use canu::Grid_Cloud;


use Carp qw(longmess cluck confess);


sub getNumberOfReadsInStore ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $nr     = 0;

    confess  if (($tag ne "all") && ($tag ne "hap") && ($tag ne "cor") && ($tag ne "obt") && ($tag ne "utg"));

    #  No file, no reads.

    return($nr)   if (! -e "./$asm.seqStore/info.txt");

    #  Read the info file.  sqStoreCreate creates this at the end.
    #
    #  The 'all' category returns the number of reads that are in the store; essentailly, the maxID of any read.
    #  The 'cor' category returns the number of reads that are available for correction.
    #  The 'obt' category returns the number of reads that are available for trimming.
    #  The 'utg' category returns the number of reads that are available for assembly.

    open(F, "< ./$asm.seqStore/info.txt") or caExit("can't open './$asm.seqStore/info.txt' for reading: $!", undef);
    while (<F>) {
        $nr = $1    if ((m/numReads\s+=\s+(\d+)/)          && ($tag eq "all"));
        $nr = $1    if ((m/numRawReads\s+=\s+(\d+)/)       && ($tag eq "cor" || $tag eq "hap"));
        $nr = $1    if ((m/numCorrectedReads\s+=\s+(\d+)/) && ($tag eq "obt"));
        $nr = $1    if ((m/numTrimmedReads\s+=\s+(\d+)/)   && ($tag eq "utg"));
    }
    close(F);

    return($nr);
}



sub getNumberOfBasesInStore ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $nb     = 0;

    confess  if (($tag ne "hap") && ($tag ne "cor") && ($tag ne "obt") && ($tag ne "utg"));

    #  No file, no bases.

    return($nb)   if (! -e "./$asm.seqStore/info.txt");

    #  Read the info file.  sqStoreCreate creates this at the end.

    open(F, "< ./$asm.seqStore/info.txt") or caExit("can't open './$asm.seqStore/info.txt' for reading: $!", undef);
    while (<F>) {
        $nb = $1    if ((m/numRawBases\s+=\s+(\d+)/)       && ($tag eq "cor" || $tag eq "hap"));
        $nb = $1    if ((m/numCorrectedBases\s+=\s+(\d+)/) && ($tag eq "obt"));
        $nb = $1    if ((m/numTrimmedBases\s+=\s+(\d+)/)   && ($tag eq "utg"));
    }
    close(F);

    return($nb);
}



sub getSizeOfSequenceStore ($) {
    my $asm    = shift @_;
    my $size   = 0;
    my $idx    = "0000";

    $size += -s "./$asm.seqStore/info";
    $size += -s "./$asm.seqStore/libraries";
    $size += -s "./$asm.seqStore/reads";

    while (-e "./$asm.seqStore/blobs.$idx") {
        $size += -s "./$asm.seqStore/blobs.$idx";
        $idx++;
    }

    return(int($size / 1024 / 1024 / 1024 + 1.5));
}



sub getExpectedCoverage ($$) {
    my $tag    = shift @_;
    my $asm    = shift @_;

    return(int(getNumberOfBasesInStore($asm, $tag) / getGlobal("genomeSize")));
}



sub createSequenceStore ($$@) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    #  If the store failed to build because of input errors and warnings, rename the store and continue.
    #  Not sure how to support this in DNANexus.

    if (-e "./$asm.seqStore.ACCEPTED") {
        rename("./$asm.seqStore.ACCEPTED",          "./$asm.seqStore");
        rename("./$asm.seqStore.BUILDING.err",      "./$asm.seqStore.err");
        return;
    }

    #  If the store failed to build and the user just reruns canu, this will be triggered.  We'll
    #  skip rebuilding the store again, and report the original error message.

    if (-e "./$asm.seqStore.BUILDING") {
        print STDERR "-- WARNING:\n";
        print STDERR "-- WARNING:  Previously failed seqStore detected.\n";
        print STDERR "-- WARNING:\n";
    }

    #  Not sure how this can occur.  Possibly the user just deleted seqStore.BUILDING and restarted?

    if ((! -e "./$asm.seqStore.BUILDING") && (-e "./$asm.seqStore.ssi")) {
        print STDERR "-- WARNING:\n";
        print STDERR "-- WARNING:  Existing sequence inputs used.\n";
        print STDERR "-- WARNING:\n";
    }

    #  Fail if there are no inputs.

    caExit("no input files specified, and store not already created, I have nothing to work on!", undef)
        if (scalar(@inputs) == 0);

    #  Convert the canu-supplied reads into correct relative paths.  This is made complicated by
    #  sqStoreCreate being run in a directory one below where we are now.

    #  At the same time, check that all files exist.

    if (! -e "./$asm.seqStore.ssi") {

        #  If any of the files are links to objects, fetch the object
        #  to local disk and update the name.
        #
        #  A similar blcok is used in SequenceStore.pm and HaplotypeReads.pm (twice).

        for (my $ff=0; $ff < scalar(@inputs); $ff++) {
            my ($inType, $inPath) = split '\0', $inputs[$ff];

            if ($inPath =~ m/dnanexus:(.*)=(.*)/) {
                my $link = $1;
                my $name = $2;

                print STDERR "-- Fetch input file './$name' from object '$link'\n";

                fetchFileFromLink($link, $name);

                $inputs[$ff] = "$inType\0./$name";
            }
        }

        #  Build a ssi file for all the raw sequence inputs.  For simplicity, we just copy in any
        #  ssi files as is.  This documents what the store was built with, etc.

        open(F, "> ./$asm.seqStore.ssi") or caExit("cant' open './$asm.seqStore.ssi' for writing: $0", undef);

        foreach my $iii (@inputs) {
            if ($iii =~ m/^-(.*)\0(.*)$/) {
                my $tech = $1;
                my $file = $2;
                my @name = split '/', $2;
                my $name = $name[scalar(@name)-1];

                $name = $1   if ($name =~ m/(.*).[xgb][z]2{0,1}$/i);
                $name = $1   if ($name =~ m/(.*).fast[aq]$/i);
                $name = $1   if ($name =~ m/(.*).f[aq]$/i);

                print F "########################################\n";
                print F "#  $tech: $file\n";
                print F "#\n";
                print F "name   $name\n";
                print F "preset $tech\n";
                print F "$file\n";
                print F "\n";

            } elsif (-e $iii) {
                print F "########################################\n";
                print F "#  $iii\n";
                print F "#\n";
                open(I, "< $iii") or caExit("can't open sqStoreCreate input '$iii' for reading: $0", undef);
                while (<I>) {
                    print F $_;
                }
                close(I);
                print F "\n";

            } else {
                caExit("unrecognized sqStoreCreate input file '$iii'", undef);
            }
        }

        close(F);
    }

    #  Load the store.

    if (! -e "./$asm.seqStore.BUILDING") {
        my $cmd;
        $cmd .= "$bin/sqStoreCreate \\\n";
        $cmd .= "  -o ./$asm.seqStore.BUILDING \\\n";
        $cmd .= "  -minlength "  . getGlobal("minReadLength")        . " \\\n";
        if (getGlobal("readSamplingCoverage") > 0) {
            $cmd .= "  -genomesize " . getGlobal("genomeSize")           . " \\\n";
            $cmd .= "  -coverage   " . getGlobal("readSamplingCoverage") . " \\\n";
            $cmd .= "  -bias       " . getGlobal("readSamplingBias")     . " \\\n";
        }
        $cmd .= "  ./$asm.seqStore.ssi \\\n";
        $cmd .= "> ./$asm.seqStore.BUILDING.err 2>&1";

        if (runCommand(".", $cmd) > 0) {
            caExit("sqStoreCreate failed", "./$asm.seqStore.BUILDING.err");
        }
    }

    #  Check for quality issues.

    if (-e "./$asm.seqStore.BUILDING.err") {
        my $nProblems = 0;

        open(F, "< ./$asm.seqStore.BUILDING.err");
        while (<F>) {
            $nProblems++   if (m/Check\syour\sreads/);
        }
        close(F);

        if ($nProblems > 0) {
            if (getGlobal("stopOnReadQuality")) {
                print STDERR "\n";
                print STDERR "Potential problems with your input reads were detected.\n";
                print STDERR "\n";
                print STDERR "Please review the logging in files:\n";
                print STDERR "  ", getcwd(), "/$asm.seqStore.BUILDING.err\n";
                print STDERR "  ", getcwd(), "/$asm.seqStore.BUILDING/errorLog\n";
                print STDERR "\n";
                print STDERR "If you wish to proceed, rename the store with the following command and restart canu.\n";
                print STDERR "\n";
                print STDERR "  mv ", getcwd(), "/$asm.seqStore.BUILDING \\\n";
                print STDERR "     ", getcwd(), "/$asm.seqStore.ACCEPTED\n";
                print STDERR "\n";
                print STDERR "Option stopOnReadQuality=false skips these checks.\n";
                exit(1);
            } else {
                print STDERR "--\n";
                print STDERR "-- WARNING:  Potential problems with your input reads were detected.\n";
                print STDERR "-- WARNING:\n";
                print STDERR "-- WARNING:  Please review the logging in files:\n";
                print STDERR "-- WARNING:    ", getcwd(), "/$asm.seqStore.BUILDING.err\n";
                print STDERR "-- WARNING:    ", getcwd(), "/$asm.seqStore.BUILDING/errorLog\n";
                print STDERR "-- \n";
                print STDERR "-- Proceeding with assembly because stopOnReadQuality=false.\n";
            }
        }
    }

    rename "./$asm.seqStore.BUILDING",             "./$asm.seqStore";
    rename "./$asm.seqStore.BUILDING.err",         "./$asm.seqStore.err";
}



sub generateReadLengthHistogram ($$) {
    my $tag    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();

    my $reads    = getNumberOfReadsInStore($asm, $tag);
    my $bases    = getNumberOfBasesInStore($asm, $tag);
    my $coverage = int(100 * $bases / getGlobal("genomeSize")) / 100;
    my $hist;

    my @rl;
    my @hi;

    #  Load the read lengths, find min and max lengths.

    open(F, "$bin/sqStoreDumpMetaData -S ./$asm.seqStore -reads |") or caExit("can't dump meta data from './$asm.seqStore': $!", undef);
    while (<F>) {
        next  if (m/readID/);             #  Skip the header.
        next  if (m/------/);

        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;
        my $l = 0;

        $l = $v[3]          if (($tag eq "cor") || ($tag eq "hap"));
        $l = $v[4]          if (($tag eq "obt"));
        $l = $v[6] - $v[5]  if (($tag eq "utg"));

        push @rl, $l   if ($l > 0);
    }
    close(F);

    @rl = sort { $a <=> $b } @rl;

    #  Generate PNG histograms if there are any reads.

    if ($reads > 0) {
        my $gnuplot = getGlobal("gnuplot");
        my $format  = getGlobal("gnuplotImageFormat");

        if ($gnuplot) {
            open(F, "> ./$asm.seqStore/readlengths-$tag.gp") or caExit("can't open './$asm.seqStore/readlengths-$tag.gp' for writing: $!", undef);
            print F "set title 'read length'\n";
            print F "set xlabel 'read length, bin width = 250'\n";
            print F "set ylabel 'number of reads'\n";
            print F "\n";
            print F "binwidth=250\n";
            print F "set boxwidth binwidth\n";
            print F "bin(x,width) = width*floor(x/width) + binwidth/2.0\n";
            print F "\n";
            print F "set terminal $format size 1024,1024\n";
            print F "set output './$asm.seqStore/readlengths-$tag.$format'\n";
            print F "plot [] './$asm.seqStore/readlengths-$tag.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
            close(F);

            open(F, "> ./$asm.seqStore/readlengths-$tag.dat") or caExit("can't open './$asm.seqStore/readlengths-$tag.dat' for writing: $!", undef);
            foreach my $rl (@rl) {
                print F "$rl\n";
            }
            close(F);

            if (runCommandSilently(".", "$gnuplot < /dev/null ./$asm.seqStore/readlengths-$tag.gp > /dev/null 2>&1", 0)) {
                print STDERR "--\n";
                print STDERR "-- WARNING: gnuplot failed.\n";
                print STDERR "--\n";
                print STDERR "----------------------------------------\n";
            }
        }
    }

    #  Generate the ASCII histogram.

    $hist  = "--\n";
    $hist .= "-- In sequence store './$asm.seqStore':\n";
    $hist .= "--   Found $reads reads.\n";
    $hist .= "--   Found $bases bases ($coverage times coverage).\n";

    if ($reads > 0) {
        my $minLen     = $rl[ 0];
        my $maxLen     = $rl[-1];
        my $scale      = 0;
        my $bucketSize = 0;

        #  Buckets of size 1000 are easy to interpret, but sometimes not ideal.

        $bucketSize = 10000;
        $bucketSize = 5000    if ($maxLen - $minLen < 1000000);
        $bucketSize = 1000    if ($maxLen - $minLen < 100000);
        $bucketSize = 100     if ($maxLen - $minLen < 10000);

        #  Generate the histogram (int truncates)

        my $mBgn = int($minLen / $bucketSize);
        my $mEnd = int($maxLen / $bucketSize);

        foreach my $rl (@rl) {
            my $b = int($rl / $bucketSize);

            $hi[$b]++;
        }

        for (my $ii=$mBgn; $ii<=$mEnd; $ii++) {                           #  Scale the *'s so that the longest has 70 of 'em
            $scale = $hi[$ii] / 70   if ($scale < $hi[$ii] / 70);
        }

        #  Draw the histogram.

        $hist .= "--\n";
        $hist .= "--   Read length histogram (one '*' equals " . int(100 * $scale) / 100 . " reads):\n";

        for (my $ii=$mBgn; $ii<=$mEnd; $ii++) {
            my $s = $ii * $bucketSize;
            my $e = $ii * $bucketSize + $bucketSize - 1;

            $hi[$ii] += 0;  #  Otherwise, cells with no count print as null.

            $hist .= sprintf("--   %6d %6d %6d %s\n", $s, $e, $hi[$ii], "*" x int($hi[$ii] / $scale));
        }
    }

    #  Abort if the read coverage is too low.

    my $minCov = getGlobal("stopOnLowCoverage");

    if ($coverage < $minCov) {
        print STDERR "--\n";
        print STDERR "-- ERROR:  Read coverage ($coverage) is too low to be useful.\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  This could be caused by an incorrect genomeSize or poor quality reads that could not\n";
        print STDERR "-- ERROR:  be sufficiently corrected.\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  You can force Canu to continue by decreasing parameter stopOnLowCoverage=$minCov,\n";
        print STDERR "-- ERROR:  however, the quality of corrected reads and/or contiguity of contigs will be poor.\n";
        print STDERR "-- \n";

        caExit("", undef);
    }

    #  Return the ASCII histogram.

    return($hist);
}



sub checkSequenceStore ($$@) {
    my $base;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    $base = "haplotype"   if ($tag eq "hap");
    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    #  Try fetching the store from object storage.  This might not be needed in all cases (e.g.,
    #  between mhap precompute and mhap compute), but it greatly simplifies stuff, like immediately
    #  here needing to check if the store exists.

    fetchSeqStore($asm);

    #  We cannot abort this step anymore.  If trimming is skipped, we need ro promote the corrected
    #  reads to trimmed reads, and also re-stash the store.
    #
    #goto allDone    if (getNumberOfReadsInStore($asm, $tag) > 0);

    #  Create the store.
    #
    #  If all goes well, we get asm.seqStore.
    #
    #  If not, we could end up with asm.seqStore.BUILDING and ask the user to examine it and rename
    #  it to asm.seqStore.ACCEPTED and restart.  On the restart, sqStoreCreateStore() detects the
    #  'ACCPETED' store and renames to asm.seqStore.

    my $histAndStash = 0;

    if (! -e "./$asm.seqStore") {
        createSequenceStore($base, $asm, @inputs);

        $histAndStash = 1;
    }

    #  Dump the list of libraries.  Various parts use this for various stuff.

    if (! -e "./$asm.seqStore/libraries.txt") {
        if (runCommandSilently(".", "$bin/sqStoreDumpMetaData -S ./$asm.seqStore -libs > ./$asm.seqStore/libraries.txt 2> /dev/null", 1)) {
            caExit("failed to generate list of libraries in store", undef);
        }
    }

    #  Query how many reads we have.

    my $nCor = getNumberOfReadsInStore($asm, "cor");   #  Number of corrected reads ready for OBT.
    my $nOBT = getNumberOfReadsInStore($asm, "obt");   #  Number of corrected reads ready for OBT.
    my $nAsm = getNumberOfReadsInStore($asm, "utg");   #  Number of trimmed reads ready for assembly.

    #  Refuse to assemble uncorrected reads.

    if (($tag eq "utg") &&    #  Assembling.
        ($nCor >  0) &&       #  And raw reads exist.
        ($nOBT == 0) &&       #  But no corrected reads.
        ($nAsm == 0)) {       #  And no trimmed reads!

        print STDERR "-- Unable to assemble uncorrected reads.\n";

        caExit("unable to assemble uncorrected reads", undef);
    }

    #  Promote corrected reads to trimmed reads, if needed.

    if (($tag eq "utg") &&    #  Assembling.
        ($nOBT > 0) &&        #  Corrected reads exist.
        ($nAsm == 0)) {       #  But no trimmed reads exist.

        print STDERR "--\n";
        print STDERR "-- WARNING:  No trimmed reads found for assembly, but untrimmed reads exist.\n";
        print STDERR "-- WARNING:  Upgrading untrimmed reads to trimmed reads for assembly.\n";

        if (runCommandSilently(".", "$bin/loadTrimmedReads -S ./$asm.seqStore > ./$asm.seqStore.upgrade.err 2>&1", 1)) {
            caExit("initializing clear ranges failed", "./$asm.seqStore.upgrade.err");
        }

        unlink "./$asm.seqStore.upgrade.err";

        $nCor = getNumberOfReadsInStore($asm, "cor");   #  Number of corrected reads ready for OBT.
        $nOBT = getNumberOfReadsInStore($asm, "obt");   #  Number of corrected reads ready for OBT.
        $nAsm = getNumberOfReadsInStore($asm, "utg");   #  Number of trimmed reads ready for assembly.

        $histAndStash = 1;
    }

    #  Make a histogram and stash the (updated) store.

    if ($histAndStash) {
        addToReport("${tag}SeqStore", generateReadLengthHistogram($tag, $asm));
        stashSeqStore($asm);
    }

    #  Refuse to conitnue if there are no reads.

    if ((($tag eq "cor") && ($nCor == 0)) ||
        (($tag eq "obt") && ($nOBT == 0)) ||
        (($tag eq "utg") && ($nAsm == 0))) {

        print STDERR "--\n";
        print STDERR "-- WARNING:\n";
        print STDERR "-- WARNING:  No raw ",       "reads detected.  Cannot proceed; empty outputs generated.\n"  if ($tag eq "cor");
        print STDERR "-- WARNING:  No corrected ", "reads detected.  Cannot proceed; empty outputs generated.\n"  if ($tag eq "obt");
        print STDERR "-- WARNING:  No trimmed ",   "reads detected.  Cannot proceed; empty outputs generated.\n"  if ($tag eq "utg");
        print STDERR "-- WARNING:\n";
        print STDERR "--\n";

        generateOutputs($asm);

        return(0);
    }

  finishStage:
    generateReport($asm);

    #  DO NOT resetIteration() here.  This function runs to completion on every
    #  restart, so resetting would obliterate the iteration count.
    #resetIteration("$tag-sqStoreCreate");

  allDone:
    stopAfter("sequenceStore");

    return(1);
}
