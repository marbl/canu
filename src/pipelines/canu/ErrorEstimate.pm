
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
 #    Sergey Koren beginning on 2016-MAY-16
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::ErrorEstimate;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(estimateCorrectedError);

use strict;
use POSIX qw(floor);

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;

#  Map subset of reads to long reads with mhap.
#  Compute resulting distribution and estimate error rate

sub estimateCorrectedError ($$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $bin     = getBinDirectory();

    $wrk = "$wrk/correction";

    my $path = "$wrk/3-estimator";

    # only run if we aren't done and were asked to
    goto allDone   if (skipStage($WRK, $asm, "errorEstimate") == 1);
    goto allDone   if (-e "$path/$asm.estimate.out");
    goto allDone   if (getGlobal("errorrate") > 0);

    print STDERR "--\n";
    print STDERR "-- ESTIMATOR (mhap)  (correction)\n";

    #  Mhap parameters - filterThreshold needs to be a string, else it is printed as 5e-06.
    #

    my ($numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer);

    $numHashes     	=  256;
    $minNumMatches 	=    4;
    $threshold     	=    0.80;
    $ordSketch     	= 1000;
    $ordSketchMer  	= getGlobal("${tag}MhapOrderedMerSize") + 2;
    my $filterThreshold = getGlobal("${tag}MhapFilterThreshold");
    my $merSize         = getGlobal("${tag}MhapMerSize");
    my $javaPath 	= getGlobal("java");

    print STDERR "--\n";
    print STDERR "-- PARAMETERS: hashes=$numHashes, minMatches=$minNumMatches, threshold=$threshold\n";
    print STDERR "--\n";

    make_path("$path");

    # subsample corrected reads, this assumes the fasta records are on a single line. We take some reads from the top and bottom of file to avoid sampling one library
    my $sampleSize = 1000;
    my $cmd = "gunzip -c $WRK/asm.correctedReads.fasta.gz |head -n $sampleSize > $path/subset.fasta";
    runCommandSilently($path, $cmd, 1);
    my $cmd = "gunzip -c $WRK/asm.correctedReads.fasta.gz |tail -n $sampleSize >> $path/subset.fasta";
    runCommandSilently($path, $cmd, 1);

    # now compute the overlaps
    $cmd  = "$javaPath -d64 -server -Xmx4g -jar $bin/mhap-" . getGlobal("${tag}MhapVersion") . ".jar ";
    $cmd .= "  --weighted -k $merSize --num-hashes $numHashes --num-min-matches $minNumMatches --ordered-sketch-size $ordSketch --ordered-kmer-size $ordSketchMer --ordered-kmer-size $ordSketchMer  --threshold $threshold --filter-threshold $filterThreshold --num-threads " . getGlobal("${tag}mhapThreads");
    $cmd .= " -s $path/subset.fasta -q $WRK/asm.correctedReads.fasta.gz 2> $path/$asm.mhap.err | $bin/errorEstimate -d 2 -m 0.985 -S - > $path/$asm.estimate.out 2> $path/$asm.estimate.err";
    runCommand($path, $cmd);

  allDone:
    return if (! -e "$path/$asm.estimate.out");

    my $errorRate = 0;
    open(L, "< $path/$asm.estimate.out") or caExit("can't open '$path/$asm.estimate.out' for reading: $!", undef);
    while (<L>) {
        $errorRate = sprintf "%.3f", ($_ / 2);
    }
    close(L);

    print STDERR "-- \n";
    if ($errorRate > 0.13) {
        print STDERR "-- Estimated error rate: " . ($errorRate*100) . "% > " . (0.13 * 100) . "% limit, capping it.\n";
        $errorRate = 0.13;
    } else {
       print STDERR "-- Estimated error rate: " . ($errorRate * 100) . "%.\n";
    }
    setErrorRate($errorRate);
    showErrorRates("--  ");
    print STDERR "-- \n";
}
