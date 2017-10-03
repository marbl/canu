
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
@EXPORT = qw(getNumberOfReadsInStore getNumberOfBasesInStore getExpectedCoverage sequenceFileExists generateReadLengthHistogram gatekeeper);

use strict;

use Cwd qw(getcwd);

use canu::Defaults;
use canu::Execution;
use canu::HTML;
use canu::Report;
use canu::Grid_Cloud;



sub getNumberOfReadsInStore ($$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $nr     = 0;

    #  No file, no reads.

    return($nr)   if (! -e "./$asm.gkpStore/info.txt");

    #  Read the info file.  gatekeeperCreate creates this at the end.

    open(F, "< ./$asm.gkpStore/info.txt") or caExit("can't open './$asm.gkpStore/info.txt' for reading: $!", undef);
    while (<F>) {
        $nr = $1    if (m/numReads\s+=\s+(\d+)/);
    }
    close(F);

    return($nr);
}



sub getNumberOfBasesInStore ($$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $nb     = 0;

    #  No file, no bases.

    return($nb)   if (! -e "./$asm.gkpStore/info.txt");

    #  Read the info file.  gatekeeperCreate creates this at the end.

    open(F, "< ./$asm.gkpStore/info.txt") or caExit("can't open './$asm.gkpStore/info.txt' for reading: $!", undef);
    while (<F>) {
        $nb = $1    if ((m/numRawBases\s+=\s+(\d+)/)       && ($base eq "correction"));
        $nb = $1    if ((m/numCorrectedBases\s+=\s+(\d+)/) && ($base eq "trimming"));
        $nb = $1    if ((m/numTrimmedBases\s+=\s+(\d+)/)   && ($base eq "unitigging"));
    }
    close(F);

    print STDERR "FOUND $nb BASES IN $base\n";

    return($nb);
}



sub getExpectedCoverage ($$) {
    my $base = shift @_;
    my $asm  = shift @_;

    return(int(getNumberOfBasesInStore($base, $asm) / getGlobal("genomeSize")));
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

    if (!-e "./$asm.gkpStore.gkp") {
        my $ff = undef;

        foreach my $iii (@inputs) {
            my ($type, $file) = split '\0', $iii;

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
        $cmd .= "  -minlength " . getGlobal("minReadLength") . " \\\n";
        $cmd .= "  -o ./$asm.gkpStore.BUILDING \\\n";
        $cmd .= "  ./$asm.gkpStore.gkp \\\n";
        $cmd .= "> ./$asm.gkpStore.BUILDING.err 2>&1";

        #  A little funny business to make gatekeeper not fail on read quality issues.
        #  A return code of 0 is total success.
        #  A return code of 1 means it found errors in the inputs, but finished.
        #  Anything larger is a crash.

        if (runCommand(".", $cmd) > 1) {
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
            print STDERR "\n";
            print STDERR "Gatekeeper detected problems in your input reads.  Please review the logging in files:\n";
            print STDERR "  ", getcwd(), "/./$asm.gkpStore.BUILDING.err\n";
            print STDERR "  ", getcwd(), "/./$asm.gkpStore.BUILDING/errorLog\n";

            if (getGlobal("stopOnReadQuality")) {
                print STDERR "If you wish to proceed, rename the store with the following commands and restart canu.\n";
                print STDERR "\n";
                print STDERR "  mv ", getcwd(), "/./$asm.gkpStore.BUILDING \\\n";
                print STDERR "     ", getcwd(), "/./$asm.gkpStore.ACCEPTED\n";
                print STDERR "\n";
                print STDERR "Or remove '", getcwd(), "/./' and re-run with stopOnReadQuality=false\n";
                print STDERR "\n";
                exit(1);
            } else {
                print STDERR "Proceeding with assembly because stopOnReadQuality=false.\n";
            }
        }
    }

    rename "./$asm.gkpStore.BUILDING",             "./$asm.gkpStore";
    rename "./$asm.gkpStore.BUILDING.err",         "./$asm.gkpStore.err";
}



sub generateReadLengthHistogram ($$) {
    my $base   = shift @_;
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

    open(F, "> ./$asm.gkpStore/readlengths-$base.dat") or caExit("can't open './$asm.gkpStore/readlengths-$base.dat' for writing: $!", undef);
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

    open(F, "> ./$asm.gkpStore/readlengths-$base.gp") or caExit("can't open './$asm.gkpStore/readlengths-$base.gp' for writing: $!", undef);
    print F "set title 'read length'\n";
    print F "set xlabel 'read length, bin width = 250'\n";
    print F "set ylabel 'number of reads'\n";
    print F "\n";
    print F "binwidth=250\n";
    print F "set boxwidth binwidth\n";
    print F "bin(x,width) = width*floor(x/width) + binwidth/2.0\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.gkpStore/readlengths-$base.lg.$format'\n";
    print F "plot [] './$asm.gkpStore/readlengths-$base.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
    print F "\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.gkpStore/readlengths-$base.sm.$format'\n";
    print F "plot [] './$asm.gkpStore/readlengths-$base.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
    close(F);

    if (runCommandSilently($base, "$gnuplot ./$asm.gkpStore/readlengths-$base.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

    #  Generate the ASCII histogram

    my $reads    = getNumberOfReadsInStore($base, $asm);
    my $bases    = getNumberOfBasesInStore($base, $asm);
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

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    #  Try fetching the store from object storage.  This might not be needed in all cases (e.g.,
    #  between mhap precompute and mhap compute), but it greatly simplifies stuff, like immediately
    #  here needing to check if the store exists.

    fetchStore("./$asm.gkpStore");

    #  If the store exists, make a symlink and be done (well, be done once we call
    #  getNumberOfReadsInStore() below).

    if ((-e "./$asm.gkpStore") && (! -e "$base/$asm.gkpStore/info")) {
        symlink("../$asm.gkpStore", "$base/$asm.gkpStore");
        return;
    }

    #  Store with reads?  Yay!  Report it, then skip.

    goto allDone    if (skipStage($asm, "$tag-gatekeeper") == 1);
    goto allDone    if (getNumberOfReadsInStore($base, $asm) > 0);

    #  Create the store.  If all goes well, we get asm.gkpStore.  If not, we could end up with
    #  asm.gkpStore.BUILDING and ask the user to examine it and rename it to asm.gkpStore.ACCEPTED
    #  and restart.  On the restart, gatekeeperCreateStore() detects the 'ACCPETED' store and
    #  renames to asm.gkpStore.

    gatekeeperCreateStore($base, $asm, @inputs)                  if (! -e "./$asm.gkpStore");

    #  Make a link to the freshly created store.

    if ((-e "./$asm.gkpStore") && (! -e "$base/$asm.gkpStore/info")) {
        symlink("../$asm.gkpStore", "$base/$asm.gkpStore");
    }

    caExit("gatekeeper store exists, but contains no reads", undef)   if (getNumberOfReadsInStore($base, $asm) == 0);

    #  Dump the list of libraries.  Various parts use this for various stuff.

    if (runCommandSilently($base, "$bin/gatekeeperDumpMetaData -G ./$asm.gkpStore -libs > ./$asm.gkpStore/libraries.txt 2> /dev/null", 1)) {
        caExit("failed to generate list of libraries in store", undef);
    }

    #  Generate a histogram of the reads.

    addToReport("${tag}GkpStore", generateReadLengthHistogram($base, $asm));

    #  Now that all the extra data is generated, stash the store.

    stashStore("./$asm.gkpStore");

  finishStage:
    emitStage($asm, "$tag-gatekeeper");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("gatekeeper");
}
