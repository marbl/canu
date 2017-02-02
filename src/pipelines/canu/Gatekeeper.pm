
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
@EXPORT = qw(getMaxReadInStore getNumberOfReadsInStore getNumberOfBasesInStore getExpectedCoverage sequenceFileExists gatekeeper);

use strict;

use canu::Defaults;
use canu::Execution;
use canu::HTML;


sub storeExists ($$) {
    my $base   = shift @_;
    my $asm    = shift @_;

    return (-e "base/$asm.gkpStore");
}


sub getMaxReadInStore ($$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $nr     = 0;

    #  No file, no reads.

    return($nr)   if (! -e "$base/$asm.gkpStore/readlengthhistogram.txt");

    #  Read the info file.  gatekeeperCreate creates this at the end.

    open(F, "< $base/$asm.gkpStore/readlengthhistogram.txt") or caExit("can't open '$base/$asm.gkpStore/readlengthhistogram.txt' for reading: $!", undef);
    while (<F>) {
       my @v = split '\s+', $_;
       $nr = $v[1];
    }
    close(F);

    return($nr);
}

sub getNumberOfReadsInStore ($$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $nr     = 0;

    #  No file, no reads.

    return($nr)   if (! -e "$base/$asm.gkpStore/info.txt");

    #  Read the info file.  gatekeeperCreate creates this at the end.

    open(F, "< $base/$asm.gkpStore/info.txt") or caExit("can't open '$base/$asm.gkpStore/info.txt' for reading: $!", undef);
    while (<F>) {
        if (m/numReads\s+=\s+(\d+)/) {
            $nr = $1;
        }
    }
    close(F);

    return($nr);
}



sub getNumberOfBasesInStore ($$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $nb     = 0;
    my $bin    = getBinDirectory();

    #  No file, no bases.

    return($nb)   if (! -e "$base/$asm.gkpStore/info.txt");

    #  Read the info file.  gatekeeperCreate creates this at the end.

    open(F, "< $base/$asm.gkpStore/reads.txt") or caExit("can't open '$base/$asm.gkpStore/reads.txt' for reading: $!", undef);
    while (<F>) {
        my @v = split '\s+', $_;
        $nb += $v[2];
    }
    close(F);

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
            return("$p$s$c")  if (-e "$p$s$c");
        }
    }

    return(undef);
}



sub gatekeeperCreateStore ($$$@) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    #  If the store failed to build because of input errors and warnings, rename the store and continue.

    if (-e "$base/$asm.gkpStore.ACCEPTED") {
        rename("$base/$asm.gkpStore.ACCEPTED",          "$base/$asm.gkpStore");
        rename("$base/$asm.gkpStore.BUILDING.err",      "$base/$asm.gkpStore.err");
        return;
    }

    #  Fail if there are no inputs.

    caExit("no input files specified, and store not already created, I have nothing to work on!", undef)
        if (scalar(@inputs) == 0);

    #  Make sure all the inputs are here.

    my $failedFiles = undef;

    foreach my $iii (@inputs) {
        my $file = $iii;  #  This stupid foreach works by reference!

        $file = $2  if ($file =~ m/^(.*)\0(.*)/);   #  Handle the raw sequence inputs.

        if (! -e $file) {
            if (defined($failedFiles)) {
                $failedFiles .= "; '$file' not found in gatekeeperCreateStore()";
            } else {
                $failedFiles = "'$file' not found in gatekeeperCreateStore()";
            }
        }
    }
    caExit($failedFiles, undef) if defined($failedFiles);

    #  Build a gkp file for all the raw sequence inputs.  For simplicity, we just copy in any gkp
    #  files as is.  This documents what gatekeeper was built with, etc.

    open(F, "> $base/$asm.gkpStore.gkp") or caExit("cant' open '$base/$asm.gkpStore.gkp' for writing: $0", undef);

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

    #  Load the store.

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

    if (runCommand($base, $cmd) > 1) {
        caExit("gatekeeper failed", "$base/$asm.gkpStore.BUILDING.err");
    }

    #  Check for quality issues.

    if (-e "$base/$asm.gkpStore.BUILDING.err") {
        my $nProblems = 0;

        open(F, "< $base/$asm.gkpStore.BUILDING.err");
        while (<F>) {
            $nProblems++   if (m/Check\syour\sreads/);
        }
        close(F);

        if ($nProblems > 0) {
            print STDERR "Gatekeeper detected problems in your input reads.  Please review the logging in files:\n";
            print STDERR "  ", getcwd(), "/$base/$asm.gkpStore.BUILDING.err\n";
            print STDERR "  ", getcwd(), "/$base/$asm.gkpStore.BUILDING/errorLog\n";

            if (getGlobal("stopOnReadQuality")) {
                print STDERR "If you wish to proceed, rename the store with the following commands and restart canu.\n";
                print STDERR "\n";
                print STDERR "  mv ", getcwd(), "/$base/$asm.gkpStore.BUILDING \\\n";
                print STDERR "     ", getcwd(), "/$base/$asm.gkpStore.ACCEPTED\n";
                print STDERR "\n";
                print STDERR "Or remove '", getcwd(), "/$base/' and re-run with stopOnReadQuality=false\n";
                print STDERR "\n";
                exit(1);
            } else {
                print STDERR "Proceeding with assembly because stopOnReadQuality=false.\n";
            }
        }
    }

    rename "$base/$asm.gkpStore.BUILDING",             "$base/$asm.gkpStore";
    rename "$base/$asm.gkpStore.BUILDING.err",         "$base/$asm.gkpStore.err";
}



sub gatekeeperGenerateReadsList ($$$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();

    if (runCommandSilently($base, "$bin/gatekeeperDumpMetaData -G ./$asm.gkpStore -reads > ./$asm.gkpStore/reads.txt 2> /dev/null", 1)) {
        caExit("failed to generate list of reads in store", undef);
    }
}

sub gatekeeperGenerateLibrariesList ($$$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();

    if (runCommandSilently($base, "$bin/gatekeeperDumpMetaData -G ./$asm.gkpStore -libs > ./$asm.gkpStore/libraries.txt 2> /dev/null", 1)) {
        caExit("failed to generate list of libraries in store", undef);
    }
}

sub gatekeeperGenerateReadLengths ($$$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();

        my $nb = 0;
        my @rl;
        my @hi;
        my $mm;

        open(F, "< $base/$asm.gkpStore/reads.txt") or caExit("can't open '$base/$asm.gkpStore/reads.txt' for reading: $!", undef);
        while (<F>) {
            my @v = split '\s+', $_;

            push @rl, $v[2];           #  Save the length
            $nb += $v[2];              #  Sum the bases
            $hi[int($v[2] / 1000)]++;  #  Add to the histogram (int truncates)
        }
        close(F);

        @rl = sort { $a <=> $b } @rl;
        $mm = int($rl[scalar(@rl)-1] / 1000);  #  max histogram value

        open(F, "> $base/$asm.gkpStore/readlengths.txt") or caExit("can't open '$base/$asm.gkpStore/readlengths.txt' for writing: $!", undef);
        foreach my $rl (@rl) {
            print F "$rl\n";
        }
        close(F);

        open(F, "> $base/$asm.gkpStore/readlengthhistogram.txt") or caExit("can't open '$base/$asm.gkpStore/readlengthhistogram.txt' for writing: $!", undef);
        for (my $ii=0; $ii<=$mm; $ii++) {
            my $s = $ii * 1000;
            my $e = $ii * 1000 + 999;

            $hi[$ii] += 0;  #  Otherwise, cells with no count print as null.

            print F "$s\t$e\t$hi[$ii]\n";
        }
        close(F);
}



sub gatekeeperGenerateReadLengthPlot ($$$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();

    my $gnuplot = getGlobal("gnuplot");
    my $format  = getGlobal("gnuplotImageFormat");

    open(F, "> $base/$asm.gkpStore/readlengths.gp") or caExit("can't open '$base/$asm.gkpStore/readlengths.gp' for writing: $!", undef);
    print F "set title 'read length'\n";
    print F "set xlabel 'read length, bin width = 250'\n";
    print F "set ylabel 'number of reads'\n";
    print F "\n";
    print F "binwidth=250\n";
    print F "set boxwidth binwidth\n";
    print F "bin(x,width) = width*floor(x/width) + binwidth/2.0\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.gkpStore/readlengths.lg.$format'\n";
    print F "plot [] './$asm.gkpStore/readlengths.txt' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
    print F "\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.gkpStore/readlengths.sm.$format'\n";
    print F "plot [] './$asm.gkpStore/readlengths.txt' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
    close(F);

    if (runCommandSilently($base, "$gnuplot ./$asm.gkpStore/readlengths.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }
}



sub gatekeeperReportReadLengthHistogram ($$$) {
    my $base   = shift @_;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();

    my $reads    = getNumberOfReadsInStore($base, $asm);
    my $bases    = getNumberOfBasesInStore($base, $asm);
    my $coverage = int(100 * $bases / getGlobal("genomeSize")) / 100;
    my $maxhist  = 0;

    open(F, "< $base/$asm.gkpStore/readlengthhistogram.txt") or caExit("can't open '$base/$asm.gkpStore/readlengthhistogram.txt' for reading: $!", undef);
    while (<F>) {
        my @v = split '\s+', $_;
        $maxhist = ($maxhist < $v[2]) ? $v[2] : $maxhist;
    }
    close(F);

    my $scale = $maxhist / 70;

    print STDERR "--\n";
    print STDERR "-- In gatekeeper store '$base/$asm.gkpStore':\n";
    print STDERR "--   Found $reads reads.\n";
    print STDERR "--   Found $bases bases ($coverage times coverage).\n";
    print STDERR "--\n";
    print STDERR "--   Read length histogram (one '*' equals ", int(100 * $scale) / 100, " reads):\n";

    open(F, "< $base/$asm.gkpStore/readlengthhistogram.txt") or caExit("can't open '$base/$asm.gkpStore/readlengthhistogram.txt' for reading: $!", undef);
    while (<F>) {
        my @v = split '\s+', $_;

        printf STDERR "--   %6d %6d %6d %s\n", $v[0], $v[1], $v[2], "*" x int($v[2] / $scale);
    }
    close(F);
}



sub gatekeeper ($$@) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    my $base;

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    #  An empty store?  Remove it and try again.

    if ((storeExists($base, $asm)) && (getNumberOfReadsInStore($base, $asm) == 0)) {
        print STDERR "-- Removing empty or incomplate gkpStore '$base/$asm.gkpStore'\n";
        runCommandSilently($base, "rm -rf ./$asm.gkpStore", 1);
    }

    #  Store with reads?  Yay!  Report it, then skip.

    goto allDone    if (skipStage($asm, "$tag-gatekeeper") == 1);
    goto allDone    if (getNumberOfReadsInStore($base, $asm) > 0);

    gatekeeperCreateStore($base, $asm, $tag, @inputs)                  if (! -e "$base/$asm.gkpStore");

    caExit("gatekeeper store exists, but contains no reads", undef)   if (getNumberOfReadsInStore($base, $asm) == 0);

    gatekeeperGenerateReadsList($base, $asm, $tag)                     if (! -e "$base/$asm.gkpStore/reads.txt");
    gatekeeperGenerateLibrariesList($base, $asm, $tag)                 if (! -e "$base/$asm.gkpStore/libraries.txt");
    gatekeeperGenerateReadLengths($base, $asm, $tag)                   if (! -e "$base/$asm.gkpStore/readlengths.txt");
    gatekeeperGenerateReadLengthPlot($base, $asm, $tag)                if (! -e "$base/$asm.gkpStore/readlengths.gp");

  finishStage:
    gatekeeperReportReadLengthHistogram($base, $asm, $tag);

    emitStage($asm, "$tag-gatekeeper");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("gatekeeper");
}
