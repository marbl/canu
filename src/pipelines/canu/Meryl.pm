
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
@EXPORT = qw(getGenomeCoverage merylConfigure merylCheck merylProcess);

use strict;

use File::Path qw(make_path);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::ErrorEstimate;
use canu::HTML;

sub getGenomeCoverage($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $merSize = shift @_;
    my $bin     = getBinDirectory();

    my $gs=`cat $wrk/0-mercounts/$asm.ms$merSize.estMerThresh.err | grep "Guessed X coverage"|awk '{print \$NF}'`;
    chomp $gs;

    # if we couldn't find the coverage, just take # bases divided by user supplied genome size as a guestimate
    if ($gs == "") {
        open(F, "$bin/gatekeeperDumpMetaData -stats -G $wrk/$asm.gkpStore | ") or caFailure("failed to read gatekeeper stats fromfrom '$wrk/$asm.gkpStore'", undef);
        while (<F>) {
           my ($junk1, $library, $junk2, $reads, $junk3, $junk4, $bases, $junk5, $average, $junk6, $min, $junk7, $max) = split '\s+', $_;
           if ($library == 0) {
              $gs = $bases / getGlobal("genomeSize");
              last;
           }
        }
       close(F);
    }

    return $gs;
}





#  Generates
#    meryl mcdat/mcidx
#    mer histogram file
#    mer histogram plots
#    mer threshold
#    frequent mers
#

#  Threshold:  Three methods to pick it.
#    Threshold  - 'auto', 'auto * X', 'auto / X', or an integer value
#    Distinct   - by the fraction distinct retained
#    Total      - by the fraction total retained
#

#  We always compute canonical mers, that are not compressed.


#  Generates   $wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta
#  stopBefore  meryl (stops before meryl itself runs)
#  stopAfter   meryl (stops after output is generated, even if it is just a symlink)


sub plotHistogram ($$$$) {
    my $wrk    = shift @_;  #  Local work directory
    my $ofile  = shift @_;
    my $suffix = shift @_;
    my $size   = shift @_;

    return  if (-e "$ofile.histogram.$suffix.png");

    open(F, "> $ofile.histogram.$suffix.gp");
    print F "\n";
    print F "unset multiplot\n";
    print F "\n";
    print F "set terminal png size $size,$size\n";
    print F "set output '$ofile.histogram.$suffix.png'\n";
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

    if (runCommandSilently("$wrk/0-mercounts", "gnuplot $ofile.histogram.$suffix.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }
}



sub merylParameters ($$$) {
    my $WRK    = shift @_;  #  Root work directory (the -d option to canu)
    ## $wrk    = $WRK;      #  Local work directory
    my $asm    = shift @_;
    my $tag    = shift @_;

    my ($wrk, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile);

    #  Find a place to run stuff.

    $wrk = "$WRK/correction"  if ($tag eq "cor");
    $wrk = "$WRK/trimming"    if ($tag eq "obt");
    $wrk = "$WRK/unitigging"  if ($tag eq "utg");

    #  Decide on which set of parameters we need to be using, and make output file names.

    if (getGlobal("${tag}Overlapper") eq "ovl") {
        $merSize      = getGlobal("${tag}OvlMerSize");
        $merThresh    = getGlobal("${tag}OvlMerThreshold");
        $merScale     = 1.0;
        $merDistinct  = getGlobal("${tag}OvlMerDistinct");
        $merTotal     = getGlobal("${tag}OvlMerTotal");

        $ffile = "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta";   #  The fasta file we should be creating (ends in FASTA).
        $ofile = "$wrk/0-mercounts/$asm.ms$merSize";                      #  The meryl database 'intermediate file'.

    } elsif (getGlobal("${tag}Overlapper") eq "mhap") {
        $merSize      = getGlobal("${tag}mhapMerSize");
        $merThresh    = undef;
        $merScale     = 1.0;
        $merDistinct  = undef;
        $merTotal     = undef;

        $ffile = "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore.gz";  #  The mhap-specific file we should be creating (ends in IGNORE).
        $ofile = "$wrk/0-mercounts/$asm.ms$merSize";                         #  The meryl database 'intermediate file'.

    } elsif (getGlobal("${tag}Overlapper") eq "minimap") {
        # do nothing
        $ffile = "$wrk/0-mercounts/$asm.ms$merSize.skip";
        make_path("$wrk/0-mercounts")  if (! -d "$wrk/0-mercounts");
        touch($ffile);
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

    #print STDERR "wrk - '$wrk'\n";
    #print STDERR "merSize - '$merSize'\n";
    #print STDERR "merThresh - '$merThresh'\n";
    #print STDERR "merScale - '$merScale'\n";
    #print STDERR "merDistinct - '$merDistinct'\n";
    #print STDERR "merTotal - '$merTotal'\n";
    #print STDERR "ffile - '$ffile'\n";
    #print STDERR "ofile - '$ofile'\n";

    return($wrk, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile);
}



sub merylConfigure ($$$) {
    my $WRK    = shift @_;  #  Root work directory (the -d option to canu)
    ## $wrk    = $WRK;      #  Local work directory
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    my ($wrk, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile) = merylParameters($WRK, $asm, $tag);

    goto allDone   if (skipStage($WRK, $asm, "$tag-merylConfigure") == 1);
    goto allDone   if  (-e "$ffile");
    goto allDone   if ((-e "$ofile.mcidx") && (-e "$ofile.mcdat"));

    make_path("$wrk/0-mercounts")  if (! -d "$wrk/0-mercounts");

    #  User supplied mers?  Just symlink to them.

    if (defined(getGlobal("${tag}OvlFrequentMers"))) {
        #my $ffile = "$wrk/0-mercounts/$asm.frequentMers.fasta";
        my $sfile = getGlobal("${tag}OvlFrequentMers");

        if (! -e $ffile) {
            caFailure("${tag}OvlFrequentMers '$sfile' not found", undef)  if (! -e $sfile);
            symlink $sfile, $ffile;
        }

        goto allDone;
    }

    #  Nope, build a script.

    my $mem = int(getGlobal("merylMemory")  * 1024 * 0.8);   #  Because meryl expects megabytes, not gigabytes.
    my $thr = getGlobal("merylThreads");

    caExit("merylMemory isn't defined?", undef)   if (!defined($mem));
    caExit("merylThreads isn't defined?", undef)  if (!defined($thr));

    open(F, "> $wrk/0-mercounts/meryl.sh") or caExit("can't open '$wrk/0-mercounts/meryl.sh' for writing: $1", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "if [ -e $wrk/$asm.tigStore/seqDB.v001.tig ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "#  Purge any previous intermediate result.  Possibly not needed, but safer.\n";
    print F "\n";
    print F "rm -f $ofile.WORKING*\n";
    print F "\n";
    print F "\$bin/meryl \\\n";
    print F "  -B -C -L 2 -v -m $merSize -threads $thr -memory $mem \\\n";
    print F "  -s $wrk/$asm.gkpStore \\\n";
    print F "  -o $ofile.WORKING \\\n";
    print F "&& \\\n";
    print F "mv $ofile.WORKING.mcdat $ofile.FINISHED.mcdat \\\n";
    print F "&& \\\n";
    print F "mv $ofile.WORKING.mcidx $ofile.FINISHED.mcidx\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

  finishStage:
    emitStage($WRK, $asm, "merylConfigure");
    buildHTML($WRK, $asm, $tag);
    stopAfter("merylConfigure");

  allDone:
}



sub merylCheck ($$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    ## $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $bin    = getBinDirectory();
    my $cmd;

    my ($wrk, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile) = merylParameters($WRK, $asm, $tag);

    #  If the frequent mer file exists, don't bother running meryl.  We don't really need the
    #  databases.

    goto allDone   if  (skipStage($WRK, $asm, "$tag-meryl") == 1);
    goto allDone   if  (-e "$ffile");
    goto allDone   if ((-e "$ofile.mcidx") && (-e "$ofile.mcdat"));

    #  If FINISHED exists, meryl finished successfully.  If WORKING, it failed during output.  If
    #  nothing, it failed before output.

    if ((! -e "$ofile.FINISHED.mcdat") ||
        (! -e "$ofile.FINISHED.mcidx")) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- meryl failed.\n";
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to generate mer counts.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- Meryl attempt $attempt begins.\n";

        emitStage($WRK, $asm, "merylCheck", $attempt);
        buildHTML($WRK, $asm, $tag);

        submitOrRunParallelJob($WRK, $asm, "meryl", "$wrk/0-mercounts", "meryl", (1));
        return;
    }

  finishStage:
    print STDERR "-- Meryl finished successfully.\n";

    rename("$ofile.FINISHED.mcdat", "$ofile.mcdat");
    rename("$ofile.FINISHED.mcidx", "$ofile.mcidx");

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "merylCheck");
    buildHTML($WRK, $asm, $tag);
    stopAfter("merylCheck");

  allDone:
}



sub merylProcess ($$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    ## $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;

    my $bin    = getBinDirectory();
    my $cmd;

    my ($wrk, $merSize, $merThresh, $merScale, $merDistinct, $merTotal, $ffile, $ofile) = merylParameters($WRK, $asm, $tag);

    goto allDone   if (skipStage($WRK, $asm, "$tag-meryl") == 1);
    goto allDone   if (-e "$ffile");

    #  A special case; if the threshold is zero, we can skip the rest.

    if ((defined($merThresh))    &&
        ($merThresh ne "auto")   &&
        ($merThresh == 0)        &&
        (!defined($merDistinct)) &&
        (!defined($merTotal))) {
        touch($ffile);
        goto allDone;
    }

    #  Dump a histogram.

    if (! -e "$ofile.histogram") {
        $cmd  = "$bin/meryl -Dh -s $ofile > $ofile.histogram 2> $ofile.histogram.info";

        if (runCommand("$wrk/0-mercounts", $cmd)) {
            rename "$ofile.histogram", "$ofile.histogram.FAILED";
            caFailure("meryl histogram failed", "$ofile.histogram.info");
        }
    }

    #  Compute a threshold, if needed

    if ($merThresh eq "auto") {
        if (! -e "$ofile.estMerThresh.out") {
            $cmd  = "$bin/estimate-mer-threshold ";
            $cmd .= " -m $ofile ";
            $cmd .= " > $ofile.estMerThresh.out ";
            $cmd .= "2> $ofile.estMerThresh.err";

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                rename "$ofile.estMerThresh.out", "$ofile.estMerThresh.out.FAILED";
                caFailure("estimate-mer-threshold failed", "$ofile.estMerThresh.err");
            }
        }

        open(F, "< $ofile.estMerThresh.out") or caFailure("failed to read estimated mer threshold from '$ofile.estMerThresh.out'", undef);
        $merThresh = <F>;
        $merThresh = int($merThresh * $merScale);
        close(F);
    }

    #  Compute a threshold based on the fraction distinct or total.

    if (defined($merDistinct) || defined($merTotal)) {
        open(F, "< $ofile.histogram") or caFailure("failed to read mer histogram from '$ofile.histogram'", undef);
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

    plotHistogram($wrk, $ofile, "lg", 1024);
    plotHistogram($wrk, $ofile, "sm", 256);

    #  Generate the frequent mers for overlapper

    if ((getGlobal("${tag}Overlapper") eq "ovl") &&
        (! -e $ffile)) {
        $cmd  = "$bin/meryl -Dt -n $merThresh -s $ofile > $ffile 2> $ffile.err";

        if (runCommand("$wrk/0-mercounts", $cmd)) {
            unlink $ffile;
            caFailure("meryl failed to dump frequent mers", "$ffile.err");
        }

        unlink "$ffile.err";
    }

    #  Generate the frequent mers for mhap
    #
    #    mer                     value           numInstances  totalKmers
    #    TTTTGTTTTTTTTTTT        0.0000044602    589           132055862
    #
    #  The fraction is just $3/$4.  I assume this is used with "--filter-threshold 0.000005".

    if ((getGlobal("${tag}Overlapper") eq "mhap") &&
        (! -e $ffile)) {

        my $totalMers = 0;
        my $maxCount  = 0;

        #  Meryl reports number of distinct canonical mers, we multiply by two to get the
        #  (approximate) number of distinct mers.  Palindromes are counted twice, oh well.

        open(F, "< $ofile.histogram.info") or die "Failed to open '$ofile.histogram.info' for reading: $!\n";
        while (<F>) {
            if (m/Found\s+(\d+)\s+mers./) {
                $totalMers = 2 * $1;
            }
            if (m/Largest\s+mercount\s+is\s+(\d+)./) {
               $maxCount = $1;
            }
        }
        close(F);
        caFailure("didn't find any mers?", "$ofile.histogram.info")  if ($totalMers == 0);

        my $filterThreshold = (getGlobal("${tag}MhapSensitivity") eq "normal") ?   getGlobal("${tag}MhapFilterThreshold") :   getGlobal("${tag}MhapFilterThreshold");  #  Also set in Meryl.pm

        my $misRate  = 0.1;
        my $minCount = defined(getGlobal("${tag}MhapFilterUnique")) ? uniqueKmerThreshold($wrk, $asm, $merSize, $misRate)+1 : int($filterThreshold * $totalMers);
        my $totalToOutput = 0;
        my $totalFiltered = 0;
        open(F, "< $ofile.histogram") or die "Failed to open '$ofile.histogram' for reading: $!\n";
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

        open(F, "$bin/meryl -Dt -n $minCount -s $ofile | ")  or die "Failed to run meryl to generate frequent mers $!\n";
        open(O, "| gzip -c > $ofile.frequentMers.ignore.gz")              or die "Failed to open '$ofile.frequentMers.ignore.gz' for writing: $!\n";
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

        if (defined(getGlobal("${tag}MhapFilterUnique"))) {
           printf STDERR "-- For %s overlapping, filtering low-occurence k-mers < %d (%.2f\%) based on estimated error of %.2f\%.\n", getGlobal("${tag}Overlapper"), $minCount, $totalFiltered, 100*estimateRawError($wrk, $asm, $tag, $merSize);
        }
        printf STDERR "-- For %s overlapping, set repeat k-mer threshold to %d.\n", getGlobal("${tag}Overlapper"), int($filterThreshold * $totalMers);
    }

    #  Report the new threshold.

    if ((getGlobal("${tag}Overlapper") eq "ovl") && ($merThresh > 0) && (getGlobal("${tag}OvlMerThreshold") ne $merThresh)) {
        print STDERR "-- Reset ${tag}OvlMerThreshold from ", getGlobal("${tag}OvlMerThreshold"), " to $merThresh.\n";
        setGlobal("${tag}OvlMerThreshold", $merThresh);
    }

  finishStage:
    if (-e "$ofile.histogram.info") {
        my $numTotal    = 0;
        my $numDistinct = 0;
        my $numUnique   = 0;
        my $largest     = 0;

        open(F, "< $ofile.histogram.info") or caFailure("can't open meryl histogram information file '$ofile.histogram.info' for reading: $!\n", undef);
        while (<F>) {
            $numTotal    = $1   if (m/Found\s(\d+)\s+mers./);
            $numDistinct = $1   if (m/Found\s(\d+)\s+distinct\smers./);
            $numUnique   = $1   if (m/Found\s(\d+)\s+unique\smers./);
            $largest     = $1   if (m/Largest\smercount\sis\s(\d+)/);
        }
        close(F);

        print STDERR "--\n";
        print STDERR "-- Found $numTotal $merSize-mers; $numDistinct distinct and $numUnique unique.  Largest count $largest.\n";

    } elsif (-z $ffile) {
        print STDERR "--\n";
        print STDERR "-- Threshold zero.  No mers will be to find.\n";

    } else {
        print STDERR "--\n";
        print STDERR "-- Using frequent mers in '", getGlobal("${tag}OvlFrequentMers"), "'\n";
    }

    unlink "$ofile.mcidx"   if (getGlobal("saveMerCounts") == 0);
    unlink "$ofile.mcdat"   if (getGlobal("saveMerCounts") == 0);

    emitStage($WRK, $asm, "$tag-meryl");
    buildHTML($WRK, $asm, $tag);
    stopAfter("meryl");

  allDone:
}
