
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
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Meryl;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(meryl);

use strict;

use canu::Defaults;
use canu::Execution;
use canu::HTML;


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
    print F "set output \"$ofile.histogram.$suffix.png\"\n";
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
    #print F "set y2tics add (\"0.6765\" 0.6765)\n";
    print F "\n";
    print F "plot [0.5:1.0] \"$ofile.histogram\" using 3:4 with lines title \"Distinct-vs-Total\"\n";
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
    #print F "set y2tics add (\"0.6765\" 0.6765)\n";
    print F "\n";
    print F "plot [0.975:1.0] \"$ofile.histogram\" using 3:4 with lines title \"Distinct-vs-Total\"\n";
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
    print F "plot [0:200] \"$ofile.histogram\" using 1:2 with lines title \"Histogram\"\n";
    close(F);

    runCommandSilently("$wrk/0-mercounts", "gnuplot $ofile.histogram.$suffix.gp > /dev/null 2>&1");
}



sub meryl ($$$) {
    my $WRK    = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk    = $WRK;      #  Local work directory
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    my ($merSize, $merThresh, $merScale, $merDistinct, $merTotal);
    my ($ffile, $ofile);

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

        $ffile = "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore";  #  The mhap-specific file we should be creating (ends in IGNORE).
        $ofile = "$wrk/0-mercounts/$asm.ms$merSize";                      #  The meryl database 'intermediate file'.

    } else {
        caFailure("unknown ${tag}Overlapper '" . getGlobal("${tag}Overlapper") . "'", undef);
    }

    #  If the frequent mer file exists, don't bother running meryl.  We don't really need the
    #  databases.

    goto allDone   if (skipStage($WRK, $asm, "$tag-meryl") == 1);
    goto allDone   if (-e "$ffile");

    print STDERR "-- MERYL (correction)\n"  if ($tag eq "cor");
    print STDERR "-- MERYL (trimming)\n"    if ($tag eq "obt");
    print STDERR "-- MERYL (assembly)\n"    if ($tag eq "utg");

    #  Make a work space.

    system("mkdir $wrk/0-mercounts")  if (! -d "$wrk/0-mercounts");

    #  User supplied mers?  Just symlink to them.

    if (defined(getGlobal("${tag}OvlFrequentMers"))) {
        my $ffile = "$wrk/0-mercounts/$asm.frequentMers.fasta";
        my $sfile = getGlobal("${tag}OvlFrequentMers");

        if (! -e $ffile) {
            caFailure("${tag}OvlFrequentMers '$sfile' not found", undef)  if (! -e $sfile);
            symlink $sfile, $ffile;
        }

        goto allDone;
    }

    #  Otherwise, run meryl, and remember the new threshold.

    #  Are we auto with modifications ("auto * X") or ("auto / X")?

    if ($merThresh =~ m/auto\s*\*\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = $1;
    }

    if ($merThresh =~ m/auto\s*\/\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = 1.0 / $1;
    }

    #  And a special case; if the threshold is zero, we can skip the rest.

    #print "thresh    $merThresh\n";
    #print "distinct  $merDistinct\n";
    #print "total     $merTotal\n";

    if ((defined($merThresh))    &&
        ($merThresh ne "auto")   &&
        ($merThresh == 0)        &&
        (!defined($merDistinct)) &&
        (!defined($merTotal))) {
        touch($ffile);
        goto allDone;
    }

    #  Build the database.

    #getAllowedResources("", "meryl");

    if (! -e "$ofile.mcdat") {
        my $mem = getGlobal("merylMemory");
        my $thr = getGlobal("merylThreads");

        $mem  = int(2 * getPhysicalMemorySize() / 3)   if (!defined($mem));
        $thr  =         getNumberOfCPUs()              if (!defined($thr));

        $mem *= 1024;  #  Because meryl expects megabytes, not gigabytes.

        $cmd  = "$bin/meryl \\\n";
        $cmd .= " -B -C -L 2 -v -m $merSize -threads $thr -memory $mem \\\n";
        $cmd .= " -s $wrk/$asm.gkpStore \\\n";
        $cmd .= " -o $ofile \\\n";
        $cmd .= "> $wrk/0-mercounts/meryl.err 2>&1";

        stopBefore("meryl", $cmd);

        if (runCommand("$wrk/0-mercounts", $cmd)) {
            caFailure("meryl failed", "$wrk/0-mercounts/meryl.err");
        }
        unlink "$wrk/0-mercounts/meryl.err";
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

        #  Meryl reports number of distinct canonical mers, we multiply by two to get the
        #  (approximate) number of distinct mers.  Palindromes are counted twice, oh well.

        open(F, "< $ofile.histogram.info") or die "Failed to open '$ofile.histogram.info' for reading: $!\n";
        while (<F>) {
            if (m/Found\s+(\d+)\s+mers./) {
                $totalMers = 2 * $1;
            }
        }
        close(F);

        caFailure("didn't find any mers?", "$ofile.histogram.info")  if ($totalMers == 0);

        my $filterThreshold = (getGlobal("${tag}MhapSensitivity") eq "normal") ?   0.000005 :   0.000005;  #  Also set in Meryl.pm
        my $minCount        = int($filterThreshold * $totalMers);

        open(F, "$bin/meryl -Dt -n $minCount -s $ofile | ")  or die "Failed to run meryl to generate frequent mers $!\n";
        open(O, "> $ofile.frequentMers.ignore.unsorted")              or die "Failed to open '$ofile.frequentMers.mhap_ignore' for writing: $!\n";

        while (!eof(F)) {
            my $h = <F>;
            my $m = <F>;  chomp $m;
            my $r = reverse $m;

            $r =~ tr/ACGTacgt/TGCAtgca/;

            if ($h =~ m/^>(\d+)/) {
                printf(O "%s\t%.16f\t$1\t$totalMers\n", $m, $1 / $totalMers);
                printf(O "%s\t%.16f\t$1\t$totalMers\n", $r, $1 / $totalMers);
            }
        }

        close(O);
        close(F);
        runCommandSilently("$wrk/0-mercounts", "cat $ofile.frequentMers.ignore.unsorted |sort -T . -rgk2 > $ofile.frequentMers.ignore"); 
    }

    #  Report the new threshold.

    if ((getGlobal("${tag}Overlapper") eq "ovl") && ($merThresh > 0) && (getGlobal("${tag}OvlMerThreshold") ne $merThresh)) {
        print STDERR "-- Reset ${tag}OvlMerThreshold from ", getGlobal("${tag}OvlMerThreshold"), " to $merThresh.\n";
        setGlobal("${tag}OvlMerThreshold", $merThresh);
    }

  finishStage:
    unlink "$ofile.mcidx"   if (getGlobal("saveMerCounts") == 0);
    unlink "$ofile.mcdat"   if (getGlobal("saveMerCounts") == 0);

    emitStage($WRK, $asm, "$tag-meryl");
    buildHTML($WRK, $asm, $tag);
    stopAfter("meryl");

  allDone:

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
}
