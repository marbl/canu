
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
 #    src/pipelines/ca3g/Meryl.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-AUG-25
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-19
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Meryl;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(merylConfigure merylCheck merylProcess);

use strict;

use File::Path 2.08 qw(make_path remove_tree);
use File::Basename;
use POSIX qw(ceil);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::ErrorEstimate;
use canu::HTML;
use canu::Report;
use canu::Grid_Cloud;



sub merylGenerateHistogram ($$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $hist;

    #  We don't know $ofile from where merylGenerateHistogram is typically called (Report.pm)
    #  and so we're forced to configure every time.

    my ($base, $path, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile) = merylParameters($asm, $tag);

    return(undef)   if (! -e "$path/$ofile.histogram");
    return(undef)   if (! -e "$path/$ofile.histogram.info");

    #  Load the statistics

    my $numTotal    = 0;
    my $numDistinct = 0;
    my $numUnique   = 0;
    my $largest     = 0;

    open(F, "< $path/$ofile.histogram.info") or caFailure("can't open meryl histogram information file '$path/$ofile.histogram.info' for reading: $!\n", undef);
    while (<F>) {
        $numTotal    = $1   if (m/Found\s(\d+)\s+mers./);
        $numDistinct = $1   if (m/Found\s(\d+)\s+distinct\smers./);
        $numUnique   = $1   if (m/Found\s(\d+)\s+unique\smers./);
        $largest     = $1   if (m/Largest\smercount\sis\s(\d+)/);
    }
    close(F);

    #  Load histogram data

    my @tc;  #  Total count
    my @fu;  #  Fraction unique
    my @ft;  #  Fraction total
    my $mc;

    open(F, "< $path/$ofile.histogram");
    while (<F>) {
        my @v = split '\s+', $_;
        $tc[$v[0]] = $v[1];
        $fu[$v[0]] = $v[2];
        $ft[$v[0]] = $v[3];
        $mc        = $v[0];  #  histogram should be sorted
    }
    close(F);

    #  Prune the high-count kmers
    #
    #  In blocks of 40, extend the histogram until the average of the next block is nearly the same
    #  as the average of this block.
if (0) {
    my $lo      = 2;
    my $hi      = 3;
    my $st      = 1;
    my $aveLast = 0;
    my $aveThis = 0;

    for (my $ii=$lo; $ii<$hi; $ii++) {
        $aveThis += $tc[$ii];
    }
    $aveThis /= ($hi - $lo);
    $aveLast  = 0;

    print STDERR "aveLast $aveLast aveThis $aveThis $lo $hi INITIAL\n";

    while (($hi < $mc) &&
           ($aveThis > 2) &&
           (($aveThis < 0.90 * $aveLast) ||
            ($aveLast < 0.90 * $aveThis))) {
        $lo += $st;
        $hi += $st;
        $st += 1;

        $aveLast = $aveThis;
        $aveThis = 0;

        for (my $ii=$lo; $ii<$hi; $ii++) {
            $aveThis += $tc[$ii];
        }
        $aveThis /= ($hi - $lo);
        print STDERR "aveLast $aveLast aveThis $aveThis $lo $hi\n";
    }

    print STDERR "aveLast $aveLast aveThis $aveThis $lo $hi FINAL\n";
}

    my @TC;
    my @FU;
    my @FT;

    my $TCmax  = 0;

    my $lo = 1;
    my $hi = 2;
    my $st = 1;

    for (my $ii=0; $ii <= 40; $ii++) {
        for (my $jj=$lo; $jj < $hi; $jj++) {
            $TC[$ii] += $tc[$jj];                                      #  Sum the counts

            $FU[$ii] = ($fu[$ii] < $FU[$ii]) ? $FU[$ii] : $fu[$jj];    #  But the fractions are already cumulative,
            $FT[$ii] = ($ft[$ii] < $FT[$ii]) ? $FT[$ii] : $ft[$jj];    #  we just need to skip zeros.
        }

        if ($ii > 0) {
            $TCmax = ($TCmax < $TC[$ii]) ? $TC[$ii] : $TCmax;
        }

        $lo  = $hi;
        $hi += $st;
        $st += 1;
    }

    my $maxY   = $lo;
    my $Xscale = $TCmax / 70;

    #  Now just draw the histogram

    $hist .= "--\n";
    $hist .= "--  $merSize-mers                                                                                           Fraction\n";
    $hist .= "--    Occurrences   NumMers                                                                         Unique Total\n";

    $lo = 1;
    $hi = 2;
    $st = 1;

    for (my $ii=0; $ii<=40; $ii++) {
        my $numXs = int($TC[$ii] / $Xscale);

        if ($numXs <= 70) {
            $hist .= sprintf("--  %6d-%6d %9d %s%s %.4f %.4f\n",
                             $lo, $hi-1, $TC[$ii],
                             "*" x      ($numXs),
                             " " x (70 - $numXs), $FU[$ii], $FT[$ii]);
        } else {
            $hist .= sprintf("--  %6d-%6d %9d %s%s %.4f %.4f\n",
                             $lo, $hi-1, $TC[$ii],
                             "*" x 67,
                             "-->", $FU[$ii], $FT[$ii]);
        }

        last   if ($hi >= $maxY);

        $lo  = $hi;
        $hi += $st;
        $st += 1;
    }

    $hist .= sprintf("--\n");
    $hist .= sprintf("-- %11d (max occurrences)\n",             $largest);
    $hist .= sprintf("-- %11d (total mers, non-unique)\n",      $numTotal    - $numUnique);
    $hist .= sprintf("-- %11d (distinct mers, non-unique)\n",   $numDistinct - $numUnique);
    $hist .= sprintf("-- %11d (unique mers)\n",                 $numUnique);

    return($hist);
}




#  Threshold:  Three methods to pick it.
#    Threshold  - 'auto', 'auto * X', 'auto / X', or an integer value
#    Distinct   - by the fraction distinct retained
#    Total      - by the fraction total retained

sub merylPlotHistogram ($$$$) {
    my $path   = shift @_;
    my $ofile  = shift @_;
    my $suffix = shift @_;
    my $size   = shift @_;  #  Size of image, not merSize!

    return  if (fileExists("$path/$ofile.histogram.$suffix.gp"));

    my $gnuplot = getGlobal("gnuplot");
    my $format  = getGlobal("gnuplotImageFormat");

    fetchFile("$path/$ofile.histogram");

    open(F, "> $path/$ofile.histogram.$suffix.gp");
    print F "\n";
    print F "unset multiplot\n";
    print F "\n";
    print F "set terminal $format size $size,$size\n";
    print F "set output '$ofile.histogram.$suffix.$format'\n";
    print F "\n";
    print F "set multiplot\n";
    print F "\n";
    print F "#  Distinct-vs-total full size plot\n";
    print F "\n";
    print F "set origin 0.0,0.0\n";
    print F "set size   1.0,1.0\n";
    print F "\n";
    print F "set xrange [0.5:1.0]\n";
    print F "set yrange [0.0:1.0]\n";
    print F "\n";
    print F "unset ytics\n";
    print F "set y2tics 0.1\n";
    #print F "set y2tics add ('0.6765' 0.6765)\n";
    print F "\n";
    print F "plot [0.5:1.0] '$ofile.histogram' using 3:4 with lines title 'Distinct-vs-Total'\n";
    print F "\n";
    print F "#  Distinct-vs-total zoom in lower left corner\n";
    print F "\n";
    print F "set origin 0.05,0.10\n";
    print F "set size   0.40,0.40\n";
    print F "\n";
    print F "set xrange [0.975:1.0]\n";
    print F "set yrange [0.4:0.80]\n";
    print F "\n";
    print F "unset ytics\n";     #  ytics on the left of the plot
    print F "set y2tics 0.1\n";  #  y2tics on the right of the plot
    #print F "set y2tics add ('0.6765' 0.6765)\n";
    print F "\n";
    print F "plot [0.975:1.0] '$ofile.histogram' using 3:4 with lines title 'Distinct-vs-Total'\n";
    print F "\n";
    print F "#  Histogram in upper left corner\n";
    print F "\n";
    print F "set origin 0.05,0.55\n";
    print F "set size   0.40,0.40\n";
    print F "\n";
    print F "set xrange [0:200]\n";
    print F "set yrange [0:30000000]\n";
    print F "\n";
    print F "unset ytics\n";      #  ytics on the left of the plot
    print F "set y2tics 10e6\n";  #  y2tics on the right of the plot
    print F "unset mytics\n";
    print F "\n";
    print F "plot [0:200] '$ofile.histogram' using 1:2 with lines title 'Histogram'\n";
    close(F);

    if (runCommandSilently($path, "$gnuplot ./$ofile.histogram.$suffix.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

    stashFile("$path/$ofile.histogram.$suffix.gp");
    stashFile("$path/$ofile.histogram.$suffix.$format");
}



sub merylParameters ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;

    my ($base, $path, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile);

    #  Find a place to run stuff.

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/0-mercounts";

    #  Decide on which set of parameters we need to be using, and make output file names.

    if (getGlobal("${tag}Overlapper") eq "ovl") {
        $merSize      = getGlobal("${tag}OvlMerSize");
        $merThresh    = getGlobal("${tag}OvlMerThreshold");
        $merScale     = 1.0;
        $merDistinct  = getGlobal("${tag}OvlMerDistinct");
        $merTotal     = getGlobal("${tag}OvlMerTotal");

        $ffile = "$asm.ms$merSize.frequentMers.fasta";   #  The fasta file we should be creating (ends in FASTA).
        $ofile = "$asm.ms$merSize";                      #  The meryl database 'intermediate file'.

    } elsif (getGlobal("${tag}Overlapper") eq "mhap") {
        $merSize      = getGlobal("${tag}mhapMerSize");
        $merThresh    = undef;
        $merScale     = 1.0;
        $merDistinct  = undef;
        $merTotal     = undef;

        $ffile = "$asm.ms$merSize.frequentMers.ignore.gz";  #  The mhap-specific file we should be creating (ends in IGNORE).
        $ofile = "$asm.ms$merSize";                         #  The meryl database 'intermediate file'.

    } elsif (getGlobal("${tag}Overlapper") eq "minimap") {
        $merSize     = 0;
        $merThresh   = 0;
        $merScale    = 1.0;
        $merDistinct = undef;
        $merTotal    = undef;

        $ffile = undef;
        $ofile = undef;

    } else {
        caFailure("unknown ${tag}Overlapper '" . getGlobal("${tag}Overlapper") . "'", undef);
    }

    #  Decode the threshold.  Auto with modifications ("auto * X") or ("auto / X")?

    if ($merThresh =~ m/auto\s*\*\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = $1;
    }

    if ($merThresh =~ m/auto\s*\/\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = 1.0 / $1;
    }

    #  Return all this goodness.

    return($base, $path, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile);
}



sub merylConfigure ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    my ($base, $path, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile) = merylParameters($asm, $tag);

    goto allDone   if (skipStage($asm, "$tag-merylConfigure") == 1);
    goto allDone   if (fileExists("$path/meryl.sh"));
    goto allDone   if (!defined($ffile));
    goto allDone   if (fileExists("$path/$ffile"));
    goto allDone   if (fileExists("$path/$ofile.mcidx") && fileExists("$path/$ofile.mcdat"));

    make_path($path)  if (! -d $path);

    #  User supplied mers?  Copy them to the proper location and exit.

    my $sfile = getGlobal("${tag}OvlFrequentMers");

    if (defined($sfile) && ! -e "$path/$ffile") {
        caFailure("${tag}OvlFrequentMers '$sfile' not found", undef)  if (! -e $sfile);
        copy($sfile, "$path/$ffile");
        stashFile("$path/$ffile");
        goto allDone;
    }

    #  No filtering?  Make an empty file and exit.

    if ((defined($merThresh))    &&
        ($merThresh ne "auto")   &&
        ($merThresh == 0)        &&
        (!defined($merDistinct)) &&
        (!defined($merTotal))) {
        touch("$path/$ffile");
        stashFile("$path/$ffile");
        goto allDone;
    }

    #  Nope, build a script for computing kmer counts.

    my $mem = int(getGlobal("merylMemory")  * 1024 * 0.8);   #  Because meryl expects megabytes, not gigabytes.
    my $thr = getGlobal("merylThreads");
    my $cov = getExpectedCoverage($base, $asm);

    caExit("merylMemory isn't defined?", undef)   if (!defined($mem));
    caExit("merylThreads isn't defined?", undef)  if (!defined($thr));

    open(F, "> $path/meryl.sh") or caExit("can't open '$path/meryl.sh' for writing: $1", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("$base/$asm.gkpStore", $path);
    print F "\n";
    print F "#  Purge any previous intermediate result.  Possibly not needed, but safer.\n";
    print F "\n";
    print F "rm -f ./$ofile.WORKING*\n";
    print F "\n";
    print F "\$bin/meryl \\\n";
    print F "  -B -C -L 2 -v -m $merSize -threads $thr -memory $mem \\\n";
    print F "  -s ../$asm.gkpStore \\\n";
    print F "  -o ./$ofile.WORKING \\\n";
    print F "&& \\\n";
    print F "mv ./$ofile.WORKING.mcdat ./$ofile.mcdat \\\n";
    print F "&& \\\n";
    print F "mv ./$ofile.WORKING.mcidx ./$ofile.mcidx\n";
    print F "\n";
    print F stashFileShellCode("$path", "$ofile.mcdat", "");
    print F "\n";
    print F stashFileShellCode("$path", "$ofile.mcidx", "");
    print F "\n";
    print F "\n";
    print F "#  Dump a histogram\n";
    print F "\n";
    print F "\$bin/meryl \\\n";
    print F "  -Dh -s ./$ofile \\\n";
    print F ">  ./$ofile.histogram.WORKING \\\n";
    print F "2> ./$ofile.histogram.info \\\n";
    print F "&& \\\n";
    print F "mv -f ./$ofile.histogram.WORKING ./$ofile.histogram\n";
    print F "\n";
    print F stashFileShellCode("$path", "$ofile.histogram", "");
    print F "\n";
    print F stashFileShellCode("$path", "$ofile.histogram.info", "");
    print F "\n";
    print F "\n";
    print F "#  Compute a nice kmer threshold.\n";
    print F "\n";
    print F "\$bin/estimate-mer-threshold \\\n";
    print F "  -h ./$ofile.histogram \\\n";
    print F "  -c $cov \\\n";
    print F ">  ./$ofile.estMerThresh.out.WORKING \\\n";
    print F "2> ./$ofile.estMerThresh.err \\\n";
    print F "&& \\\n";
    print F "mv ./$ofile.estMerThresh.out.WORKING ./$ofile.estMerThresh.out\n";
    print F "\n";
    print F stashFileShellCode("$path", "$ofile.estMerThresh.out", "");
    print F "\n";
    print F stashFileShellCode("$path", "$ofile.estMerThresh.err", "");
    print F "\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

    makeExecutable("$path/meryl.sh");
    stashFile("$path/meryl.sh");

  finishStage:
    emitStage($asm, "merylConfigure");
    buildHTML($asm, $tag);

  allDone:
}



sub merylCheck ($$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $bin     = getBinDirectory();
    my $cmd;

    my ($base, $path, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile) = merylParameters($asm, $tag);

    #  If the frequent mer file exists, don't bother running meryl.  We don't really need the
    #  databases.

    goto allDone      if (skipStage($asm, "$tag-meryl") == 1);
    goto allDone      if (fileExists("$path/meryl.success"));
    goto finishStage  if (!defined($ffile));
    goto finishStage  if (fileExists("$path/$ffile"));
    goto finishStage  if (fileExists("$path/$ofile.mcidx") && fileExists("$path/$ofile.mcdat"));

    fetchFile("$path/meryl.sh");

    #  Since there is only one job, if we get here, we're not done.  Any other 'check' function
    #  shows how to process multiple jobs.  This only checks for the existence of the final outputs.
    #  (unitigger is the same)

    #  If too many attempts, give up.

    if ($attempt >= getGlobal("canuIterationMax")) {
        print STDERR "--\n";
        print STDERR "-- Meryl failed, tried $attempt times, giving up.\n";
        print STDERR "--\n";
        caExit(undef, undef);
    }

    if ($attempt > 0) {
        print STDERR "--\n";
        print STDERR "-- Meryl failed, retry.\n";
        print STDERR "--\n";
    }

    #  Otherwise, run some jobs.

    emitStage($asm, "merylCheck", $attempt);
    buildHTML($asm, $tag);

    submitOrRunParallelJob($asm, "meryl", $path, "meryl", (1));
    return;

  finishStage:
    print STDERR "-- Meryl finished successfully.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/meryl.success") or caExit("can't open '$path/meryl.success' for writing: $!", undef);
    close(F);

    stashFile("$path/meryl.success");

    emitStage($asm, "merylCheck");
    buildHTML($asm, $tag);

  allDone:
}



sub merylProcess ($$) {
    my $asm     = shift @_;
    my $tag     = shift @_;

    my $bin     = getBinDirectory();
    my $cmd;

    my ($base, $path, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile) = merylParameters($asm, $tag);

    #  ffile exists if we've already output it here, or if user supplied a file, or if user wants no masking.

    goto allDone   if (skipStage($asm, "$tag-meryl") == 1);
    goto allDone   if (fileExists("$path/$ffile"));

    #  Compute a threshold, if needed.

    if ($merThresh eq "auto") {
        fetchFile("$path/$ofile.estMerThresh.out");

        open(F, "< $path/$ofile.estMerThresh.out") or caFailure("failed to read estimated mer threshold from '$path/$ofile.estMerThresh.out'", undef);
        $merThresh = <F>;
        $merThresh = int($merThresh * $merScale) + 1;
        close(F);
    }

    #  Compute a threshold based on the fraction distinct or total.

    if (defined($merDistinct) || defined($merTotal)) {
        fetchFile("$path/$ofile.histogram");

        open(F, "< $path/$ofile.histogram") or caFailure("failed to read mer histogram from '$path/$ofile.histogram'", undef);
        while (<F>) {
            my ($threshold, $num, $distinct, $total) = split '\s+', $_;

            if (($merThresh > 0) && ($merThresh < $threshold)) {
                print STDERR "-- Supplied merThreshold $merThresh is the smallest.\n";
                last;
            }

            if ((defined($merDistinct)) && ($merDistinct <= $distinct)) {
                $merThresh = (($merThresh > 0) && ($merThresh < $threshold)) ? $merThresh : $threshold;
                print STDERR "-- Supplied merDistinct $merDistinct with threshold $threshold is the smallest.\n";
                last;
            }

            if ((defined($merTotal)) && ($merTotal <= $total)) {
                $merThresh = (($merThresh > 0) && ($merThresh < $threshold)) ? $merThresh : $threshold;
                print STDERR "-- Supplied merTotal $merTotal with threshold $threshold is the smallest.\n";
                last;
            }
        }
        close(F);
    }

    #  Plot the histogram - annotated with the thesholds

    merylPlotHistogram($path, $ofile, "lg", 1024);    #  $ofile has merSize encoded in it
    merylPlotHistogram($path, $ofile, "sm", 256);

    #  Display the histogram, and save to the report.  Shouldn't this (and the plots above)
    #  go in finishStage?

    addToReport("${tag}Meryl", merylGenerateHistogram($asm, $tag));

    #  Generate the frequent mers for overlapper

    if (getGlobal("${tag}Overlapper") eq "ovl") {
        fetchFile("$path/$ofile.mcdat");
        fetchFile("$path/$ofile.mcidx");

        if ((! -e "$path/$ofile.mcdat") ||
            (! -e "$path/$ofile.mcdat")) {
            caFailure("meryl can't dump frequent mers, databases don't exist.  Remove $path/meryl.success to try again.", undef);
        }

        if (runCommand($path, "$bin/meryl -Dt -n $merThresh -s ./$ofile > ./$ffile 2> ./$ffile.err")) {
            unlink "$path/$ffile";
            caFailure("meryl failed to dump frequent mers", "$path/$ffile.err");
        }
        unlink "$path/$ffile.err";

        stashFile("$path/$ffile");
    }

    #  Generate the frequent mers for mhap
    #
    #    mer                     value           numInstances  totalKmers
    #    TTTTGTTTTTTTTTTT        0.0000044602    589           132055862
    #
    #  The fraction is just $3/$4.  I assume this is used with "--filter-threshold 0.000005".

    if (getGlobal("${tag}Overlapper") eq "mhap") {
        my $totalMers = 0;
        my $maxCount  = 0;

        fetchFile("$path/$ofile.histogram");
        fetchFile("$path/$ofile.histogram.info");

        #  Meryl reports number of distinct canonical mers, we multiply by two to get the
        #  (approximate) number of distinct mers.  Palindromes are counted twice, oh well.

        open(F, "< $path/$ofile.histogram.info") or die "Failed to open '$path/$ofile.histogram.info' for reading: $!\n";
        while (<F>) {
            if (m/Found\s+(\d+)\s+mers./) {
                $totalMers = 2 * $1;
            }
            if (m/Largest\s+mercount\s+is\s+(\d+)./) {
               $maxCount = $1;
            }
        }
        close(F);
        caFailure("didn't find any mers?", "$path/$ofile.histogram.info")  if ($totalMers == 0);

        my $filterThreshold = getGlobal("${tag}MhapFilterThreshold");
        my $misRate         = 0.1;
        my $minCount        = int($filterThreshold * $totalMers);
        my $totalToOutput   = 0;
        my $totalFiltered   = 0;

        if (defined(getGlobal("${tag}MhapFilterUnique"))) {
            $minCount = uniqueKmerThreshold($base, $asm, $merSize, $misRate) + 1;
        }

        open(F, "< $path/$ofile.histogram") or die "Failed to open '$path/$ofile.histogram' for reading: $!\n";
        while (<F>) {
           my ($kCount, $occurences, $cumsum, $faction) = split '\s+', $_;
           if ($kCount < $minCount) {
              $totalFiltered = $cumsum * 100;
           }
           if ($kCount >= $minCount) {
              $totalToOutput += $occurences;
           }
        }
        close(F);
        $totalToOutput *= 2; # for the reverse complement

        fetchFile("$path/$ofile.mcdat");
        fetchFile("$path/$ofile.mcidx");

        open(F, "$bin/meryl -Dt -n $minCount -s $path/$ofile | ")    or die "Failed to run meryl to generate frequent mers $!\n";
        open(O, "| gzip -c > $path/$ofile.frequentMers.ignore.gz")   or die "Failed to open '$path/$ofile.frequentMers.ignore.gz' for writing: $!\n";

        printf(O "%d\n", $totalToOutput);

        while (!eof(F)) {
            my $h = <F>;
            my $m = <F>;  chomp $m;
            my $r = reverse $m;

            $r =~ tr/ACGTacgt/TGCAtgca/;

            if ($h =~ m/^>(\d+)/) {
                printf(O "%s\t%e\n", $m, $1 / $totalMers);
                printf(O "%s\t%e\n", $r, $1 / $totalMers);
            }
        }
        close(O);
        close(F);

        stashFile("$path/$ffile");

        if (defined(getGlobal("${tag}MhapFilterUnique"))) {
           printf STDERR "-- For %s overlapping, filtering low-occurence k-mers < %d (%.2f\%) based on estimated error of %.2f\%.\n", getGlobal("${tag}Overlapper"), $minCount, $totalFiltered, 100*estimateRawError($base, $asm, $tag, $merSize);
        }
        printf STDERR "-- For %s overlapping, set repeat k-mer threshold to %d.\n", getGlobal("${tag}Overlapper"), int($filterThreshold * $totalMers);
    }

    #  Report the new threshold.

    if ((getGlobal("${tag}Overlapper") eq "ovl") && ($merThresh > 0) && (getGlobal("${tag}OvlMerThreshold") ne $merThresh)) {
        print STDERR "-- Reset ${tag}OvlMerThreshold from ", getGlobal("${tag}OvlMerThreshold"), " to $merThresh.\n";
        setGlobal("${tag}OvlMerThreshold", $merThresh);
    }

  finishStage:
    fetchFile("$path/$ofile.histogram.info");
    fetchFile("$path/$ffile");

    if (-e "$path/$ofile.histogram.info") {
        my $numTotal    = 0;
        my $numDistinct = 0;
        my $numUnique   = 0;
        my $largest     = 0;

        open(F, "< $path/$ofile.histogram.info") or caFailure("can't open meryl histogram information file '$path/$ofile.histogram.info' for reading: $!\n", undef);
        while (<F>) {
            $numTotal    = $1   if (m/Found\s(\d+)\s+mers./);
            $numDistinct = $1   if (m/Found\s(\d+)\s+distinct\smers./);
            $numUnique   = $1   if (m/Found\s(\d+)\s+unique\smers./);
            $largest     = $1   if (m/Largest\smercount\sis\s(\d+)/);
        }
        close(F);

        print STDERR "--\n";
        print STDERR "-- Found $numTotal $merSize-mers; $numDistinct distinct and $numUnique unique.  Largest count $largest.\n";

    } elsif (-z "$path/$ffile") {
        print STDERR "--\n";
        print STDERR "-- Threshold zero.  No mers will be masked.\n";

    } else {
        print STDERR "--\n";
        print STDERR "-- Using frequent mers in '", getGlobal("${tag}OvlFrequentMers"), "'\n";
    }

    unlink "$path/$ofile.mcidx"   if (getGlobal("saveMerCounts") == 0);
    unlink "$path/$ofile.mcdat"   if (getGlobal("saveMerCounts") == 0);

    emitStage($asm, "$tag-meryl");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("meryl");
}
