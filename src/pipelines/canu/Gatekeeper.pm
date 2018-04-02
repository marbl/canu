
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
 #  This file is derived from:
 #
 #    src/pipelines/ca3g/Gatekeeper.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-27
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2016-FEB-29
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Gatekeeper;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getNumberOfReadsEarliestVersion
             getNumberOfReadsInStore
             getNumberOfBasesInStore
             getSizeOfGatekeeperStore
             getExpectedCoverage
             getExpectedCoverageEarliestVersion
             sequenceFileExists
             generateReadLengthHistogram
             gatekeeper);

use strict;

use Cwd qw(getcwd);

use canu::Defaults;
use canu::Execution;
use canu::Report;
use canu::Grid_Cloud;


use Carp qw(longmess cluck confess);

sub getNumberOfReadsEarliestVersion($) {
    my $asm    = shift @_;

    # get the earliest count we have in the store
    my $maxID    = getNumberOfReadsInStore("cor", $asm);
    $maxID       = getNumberOfReadsInStore("obt", $asm) if $maxID == 0;
    $maxID       = getNumberOfReadsInStore("utg", $asm) if $maxID == 0;

    return $maxID;
}

sub getNumberOfReadsInStore ($$) {
    my $tag    = shift @_;
    my $asm    = shift @_;
    my $nr     = 0;

    confess  if (($tag ne "hap") && ($tag ne "cor") && ($tag ne "obt") && ($tag ne "utg"));

    #  No file, no reads.

    return($nr)   if (! -e "./$asm.gkpStore/info.txt");

    #  Read the info file.  gatekeeperCreate creates this at the end.

    open(F, "< ./$asm.gkpStore/info.txt") or caExit("can't open './$asm.gkpStore/info.txt' for reading: $!", undef);
    while (<F>) {
        $nr = $1    if ((m/numRawReads\s+=\s+(\d+)/)       && ($tag eq "cor" || $tag eq "hap"));
        $nr = $1    if ((m/numCorrectedReads\s+=\s+(\d+)/) && ($tag eq "obt"));
        $nr = $1    if ((m/numTrimmedReads\s+=\s+(\d+)/)   && ($tag eq "utg"));
    }
    close(F);

    return($nr);
}



sub getNumberOfBasesInStore ($$) {
    my $tag    = shift @_;
    my $asm    = shift @_;
    my $nb     = 0;

    confess  if (($tag ne "hap") && ($tag ne "cor") && ($tag ne "obt") && ($tag ne "utg"));

    #  No file, no bases.

    return($nb)   if (! -e "./$asm.gkpStore/info.txt");

    #  Read the info file.  gatekeeperCreate creates this at the end.

    open(F, "< ./$asm.gkpStore/info.txt") or caExit("can't open './$asm.gkpStore/info.txt' for reading: $!", undef);
    while (<F>) {
        $nb = $1    if ((m/numRawBases\s+=\s+(\d+)/)       && ($tag eq "cor" || $tag eq "hap"));
        $nb = $1    if ((m/numCorrectedBases\s+=\s+(\d+)/) && ($tag eq "obt"));
        $nb = $1    if ((m/numTrimmedBases\s+=\s+(\d+)/)   && ($tag eq "utg"));
    }
    close(F);

    return($nb);
}



sub getSizeOfGatekeeperStore ($) {
    my $asm    = shift @_;
    my $size   = 0;
    my $idx    = "0000";

    $size += -s "./$asm.gkpStore/info";
    $size += -s "./$asm.gkpStore/libraries";
    $size += -s "./$asm.gkpStore/reads";

    while (-e "./$asm.gkpStore/blobs.$idx") {
        $size += -s "./$asm.gkpStore/blobs.$idx";
        $idx++;
    }

    return(int($size / 1024 / 1024 / 1024 + 1.5));
}



sub getExpectedCoverageEarliestVersion($) {
    my $asm    = shift @_;

    # get earliest count we can
    my $coverage     = getExpectedCoverage("cor", $asm);
    $coverage = getExpectedCoverage("obt", $asm) if $coverage == 0;
    $coverage = getExpectedCoverage("utg", $asm) if $coverage == 0;

    return $coverage;
}

sub getExpectedCoverage ($$) {
    my $tag    = shift @_;
    my $asm    = shift @_;

    return(int(getNumberOfBasesInStore($tag, $asm) / getGlobal("genomeSize")));
}



#  Returns undef if a sequence file with the supplied name cannot be found.  Common suffices and compressions are tried.
#  Otherwise, returns the found sequence file.

sub sequenceFileExists ($) {
    my $p = shift @_;

    foreach my $s ("", ".fasta", ".fastq", ".fa", ".fq") {
        foreach my $c ("", ".gz", ".xz") {
            return("$p$s$c")  if (fileExists("$p$s$c"));
        }
    }

    return(undef);
}



sub gatekeeperCreateStore ($$@) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    #  If the store failed to build because of input errors and warnings, rename the store and continue.
    #  Not sure how to support this in DNANexus.

    if (-e "./$asm.gkpStore.ACCEPTED") {
        rename("./$asm.gkpStore.ACCEPTED",          "./$asm.gkpStore");
        rename("./$asm.gkpStore.BUILDING.err",      "./$asm.gkpStore.err");
        return;
    }

    #  If the store failed to build and the user just reruns canu, this will be triggered.  We'll
    #  skip rebuilding the store again, and report the original error message.

    if (-e "./$asm.gkpStore.BUILDING") {
        print STDERR "-- WARNING:\n";
        print STDERR "-- WARNING:  Previously failed gkpStore detected.\n";
        print STDERR "-- WARNING:\n";
    }

    #  Not sure how this can occur.  Possibly the user just deleted gkpStore.BUILDING and restarted?

    if ((! -e "./$asm.gkpStore.BUILDING") && (-e "./$asm.gkpStore.gkp")) {
        print STDERR "-- WARNING:\n";
        print STDERR "-- WARNING:  Existing sequence inputs used.\n";
        print STDERR "-- WARNING:\n";
    }

    #  Fail if there are no inputs.

    caExit("no input files specified, and store not already created, I have nothing to work on!", undef)
        if (scalar(@inputs) == 0);

    #  Convert the canu-supplied reads into correct relative paths.  This is made complicated by
    #  gatekeeperCreate being run in a directory one below where we are now.

    #  At the same time, check that all files exist.

    if (! -e "./$asm.gkpStore.gkp") {
        my $ff = undef;

        foreach my $iii (@inputs) {
            my ($type, $file) = split '\0', $iii;

            #  REMOVE THIS!
            if (($file =~ m/\.correctedReads\./) ||
                ($file =~ m/\.trimmedReads\./)) {
                fetchFile($file);

                chdir($base);                                 #  Move to where we run the command
                $file = "../$file"      if (-e "../$file");   #  If file exists up one dir, it's our file
                $iii = "$type\0$file";                        #  Rewrite the option
                chdir("..");                                  #  ($file is used below too)
            }

            chdir($base);
            $ff .= (defined($ff) ? "\n  " : "") . "reads '$file' not found."  if (! -e $file);
            chdir("..");
        }

        caExit($ff, undef) if defined($ff);


        #  Build a gkp file for all the raw sequence inputs.  For simplicity, we just copy in any gkp
        #  files as is.  This documents what gatekeeper was built with, etc.

        open(F, "> ./$asm.gkpStore.gkp") or caExit("cant' open './$asm.gkpStore.gkp' for writing: $0", undef);

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
                open(I, "< $iii") or caExit("can't open gatekeeper input '$iii' for reading: $0", undef);
                while (<I>) {
                    print F $_;
                }
                close(I);
                print F "\n";

            } else {
                caExit("unrecognized gatekeeper input file '$iii'", undef);
            }
        }

        close(F);
    }

    #  Load the store.

    if (! -e "./$asm.gkpStore.BUILDING") {
        my $cmd;
        $cmd .= "$bin/gatekeeperCreate \\\n";
        $cmd .= "  -o ./$asm.gkpStore.BUILDING \\\n";
        $cmd .= "  -minlength "  . getGlobal("minReadLength")        . " \\\n";
        if (getGlobal("readSamplingCoverage") > 0) {
            $cmd .= "  -genomesize " . getGlobal("genomeSize")           . " \\\n";
            $cmd .= "  -coverage   " . getGlobal("readSamplingCoverage") . " \\\n";
            $cmd .= "  -bias       " . getGlobal("readSamplingBias")     . " \\\n";
        }
        $cmd .= "  ./$asm.gkpStore.gkp \\\n";
        $cmd .= "> ./$asm.gkpStore.BUILDING.err 2>&1";

        if (runCommand(".", $cmd) > 0) {
            caExit("gatekeeper failed", "./$asm.gkpStore.BUILDING.err");
        }
    }

    #  Check for quality issues.

    if (-e "./$asm.gkpStore.BUILDING.err") {
        my $nProblems = 0;

        open(F, "< ./$asm.gkpStore.BUILDING.err");
        while (<F>) {
            $nProblems++   if (m/Check\syour\sreads/);
        }
        close(F);

        if ($nProblems > 0) {
            if (getGlobal("stopOnReadQuality")) {
                print STDERR "\n";
                print STDERR "Gatekeeper detected potential problems in your input reads.\n";
                print STDERR "\n";
                print STDERR "Please review the logging in files:\n";
                print STDERR "  ", getcwd(), "/$asm.gkpStore.BUILDING.err\n";
                print STDERR "  ", getcwd(), "/$asm.gkpStore.BUILDING/errorLog\n";
                print STDERR "\n";
                print STDERR "If you wish to proceed, rename the store with the following command and restart canu.\n";
                print STDERR "\n";
                print STDERR "  mv ", getcwd(), "/$asm.gkpStore.BUILDING \\\n";
                print STDERR "     ", getcwd(), "/$asm.gkpStore.ACCEPTED\n";
                print STDERR "\n";
                print STDERR "Option stopOnReadQuality=false skips these checks.\n";
                exit(1);
            } else {
                print STDERR "--\n";
                print STDERR "-- WARNING:  Gatekeeper detected potential problems in your input reads.\n";
                print STDERR "-- WARNING:\n";
                print STDERR "-- WARNING:  Please review the logging in files:\n";
                print STDERR "-- WARNING:    ", getcwd(), "/$asm.gkpStore.BUILDING.err\n";
                print STDERR "-- WARNING:    ", getcwd(), "/$asm.gkpStore.BUILDING/errorLog\n";
                print STDERR "-- \n";
                print STDERR "-- Proceeding with assembly because stopOnReadQuality=false.\n";
            }
        }
    }

    rename "./$asm.gkpStore.BUILDING",             "./$asm.gkpStore";
    rename "./$asm.gkpStore.BUILDING.err",         "./$asm.gkpStore.err";
}



sub generateReadLengthHistogram ($$) {
    my $tag    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();

    my $nb = 0;
    my @rl;
    my @hi;
    my $mm;
    my $minLen = 999999;
    my $maxLen = 0;

    open(F, "$bin/gatekeeperDumpMetaData -G ./$asm.gkpStore -reads |") or caExit("can't dump meta data from './$asm.gkpStore': $!", undef);
    while (<F>) {
        next  if (m/readID/);             #  Skip the header.
        next  if (m/------/);

        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        if ($v[2] > 0) {
            push @rl, $v[2];              #  Save the length
            $nb += $v[2];                 #  Sum the bases

            $minLen = ($minLen < $v[2]) ? $minLen : $v[2];
            $maxLen = ($v[2] < $maxLen) ? $maxLen : $v[2];
        }
    }
    close(F);

    @rl = sort { $a <=> $b } @rl;

    #  Buckets of size 1000 are easy to interpret, but sometimes not ideal.

    my $bucketSize = 0;

    if      ($maxLen - $minLen < 10000) {
        $bucketSize = 100;
    } elsif ($maxLen - $minLen < 100000) {
        $bucketSize = 1000;
    } elsif ($maxLen - $minLen < 1000000) {
        $bucketSize = 5000;
    } else {
        $bucketSize = 10000;
    }

    #  Generate the histogram (int truncates)

    foreach my $rl (@rl) {
        my $b = int($rl / $bucketSize);

        $hi[$b]++;
    }

    $mm = int($maxLen / $bucketSize);  #  Max histogram value

    #  Write the sorted read lengths (for gnuplot) and the maximum read length (for correction consensus)

    open(F, "> ./$asm.gkpStore/readlengths-$tag.dat") or caExit("can't open './$asm.gkpStore/readlengths-$tag.dat' for writing: $!", undef);
    foreach my $rl (@rl) {
        print F "$rl\n";
    }
    close(F);

    #open(F, "> ./$asm.gkpStore/maxreadlength.txt") or caExit("can't open './$asm.gkpStore/maxreadlength.txt' for writing: $!", undef);
    #print F "$maxLen\n";
    #close(F);

    #  Generate PNG histograms

    my $gnuplot = getGlobal("gnuplot");
    my $format  = getGlobal("gnuplotImageFormat");

    open(F, "> ./$asm.gkpStore/readlengths-$tag.gp") or caExit("can't open './$asm.gkpStore/readlengths-$tag.gp' for writing: $!", undef);
    print F "set title 'read length'\n";
    print F "set xlabel 'read length, bin width = 250'\n";
    print F "set ylabel 'number of reads'\n";
    print F "\n";
    print F "binwidth=250\n";
    print F "set boxwidth binwidth\n";
    print F "bin(x,width) = width*floor(x/width) + binwidth/2.0\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.gkpStore/readlengths-$tag.lg.$format'\n";
    print F "plot [] './$asm.gkpStore/readlengths-$tag.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
    print F "\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.gkpStore/readlengths-$tag.sm.$format'\n";
    print F "plot [] './$asm.gkpStore/readlengths-$tag.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
    close(F);

    if (runCommandSilently(".", "$gnuplot ./$asm.gkpStore/readlengths-$tag.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

    #  Generate the ASCII histogram

    my $reads    = getNumberOfReadsInStore($tag, $asm);
    my $bases    = getNumberOfBasesInStore($tag, $asm);
    my $coverage = int(100 * $bases / getGlobal("genomeSize")) / 100;
    my $scale    = 0;
    my $hist;

    for (my $ii=0; $ii<=$mm; $ii++) {                           #  Scale the *'s so that the longest has 70 of 'em
        $scale = $hi[$ii] / 70   if ($scale < $hi[$ii] / 70);
    }

    $hist  = "--\n";
    $hist .= "-- In gatekeeper store './$asm.gkpStore':\n";
    $hist .= "--   Found $reads reads.\n";
    $hist .= "--   Found $bases bases ($coverage times coverage).\n";
    $hist .= "--\n";
    $hist .= "--   Read length histogram (one '*' equals " . int(100 * $scale) / 100 . " reads):\n";

    for (my $ii=0; $ii<=$mm; $ii++) {
        my $s = $ii * $bucketSize;
        my $e = $ii * $bucketSize + $bucketSize - 1;

        $hi[$ii] += 0;  #  Otherwise, cells with no count print as null.

        $hist .= sprintf("--   %6d %6d %6d %s\n", $s, $e, $hi[$ii], "*" x int($hi[$ii] / $scale));
    }

    return($hist);
}



sub gatekeeper ($$@) {
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

    fetchStore("./$asm.gkpStore");

    #  We cannot abort this step anymore.  If trimming is skipped, we need ro promote the corrected
    #  reads to trimmed reads, and also re-stash the store.
    #
    #goto allDone    if (skipStage($asm, "$tag-gatekeeper") == 1);
    #goto allDone    if (getNumberOfReadsInStore($tag, $asm) > 0);

    #  Create the store.
    #
    #  If all goes well, we get asm.gkpStore.
    #
    #  If not, we could end up with asm.gkpStore.BUILDING and ask the user to examine it and rename
    #  it to asm.gkpStore.ACCEPTED and restart.  On the restart, gatekeeperCreateStore() detects the
    #  'ACCPETED' store and renames to asm.gkpStore.

    my $histAndStash = 0;

    if (! -e "./$asm.gkpStore") {
        gatekeeperCreateStore($base, $asm, @inputs);

        $histAndStash = 1;
    }

    #  Dump the list of libraries.  Various parts use this for various stuff.

    if (! -e "./$asm.gkpStore/libraries.txt") {
        if (runCommandSilently(".", "$bin/gatekeeperDumpMetaData -G ./$asm.gkpStore -libs > ./$asm.gkpStore/libraries.txt 2> /dev/null", 1)) {
            caExit("failed to generate list of libraries in store", undef);
        }
    }

    #  Most of the pipeline still expects a gkpStore to exist in the stage subdirectories.  So make it exist.

    symlink("../$asm.gkpStore", "$base/$asm.gkpStore")    if ((-e "./$asm.gkpStore") && (! -e "$base/$asm.gkpStore/info"));

    if (! -e "$base/$asm.gkpStore/info") {
        print STDERR "ERROR:\n";
        print STDERR "ERROR:  Failed to create a symlink to '$asm.gkpStore' from within the '$base' directory.\n";
        print STDERR "ERROR:  This is known to happen on VirtualBox when running in the 'shared' directory.\n";
        print STDERR "ERROR:\n";

        caExit("failed to make symlink '$base/$asm.gkpStore' to '../$asm.gkpStore'", undef);
    }

    #  Query how many reads we have.

    my $nCor = getNumberOfReadsInStore("cor", $asm);   #  Number of corrected reads ready for OBT.
    my $nOBT = getNumberOfReadsInStore("obt", $asm);   #  Number of corrected reads ready for OBT.
    my $nAsm = getNumberOfReadsInStore("utg", $asm);   #  Number of trimmed reads ready for assembly.

    #  Refuse to assemble uncorrected reads (but keep going if there are no reads at all).

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

        if (runCommandSilently(".", "$bin/loadTrimmedReads -G ./$asm.gkpStore > ./$asm.gkpStore.upgrade.err 2>&1", 1)) {
            caExit("initializing clear ranges failed", "./$asm.gkpStore.upgrade.err");
        }

        unlink "./$asm.gkpStore.upgrade.err";

        $histAndStash = 1;
    }

    #  Make a histogram and stash the (updated) store.

    if ($histAndStash) {
        addToReport("${tag}GkpStore", generateReadLengthHistogram($tag, $asm));
        stashStore("./$asm.gkpStore");
    }

  finishStage:
    ;  #  Perl 5.10 (at least) is VERY unhappy about having two adjacent labels.
    #  DO NOT emitStage() here.  It resets canuIteration, and doesn't need to.
    #emitStage($asm, "$tag-gatekeeper");

  allDone:
    stopAfter("gatekeeper");
}
