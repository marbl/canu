
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
@EXPORT = qw(estimateKmerError estimateRawError estimateCorrectedError uniqueKmerThreshold);

use strict;
use POSIX qw(floor);

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;

sub fac($) {
    my $x = shift @_;

    return 1 if($x == 0);
    return 1 if($x == 1);
    return $x * fac($x - 1);
}

sub poisson_pdf ($$) {
    my $lambda = shift @_;
    my $k = shift @_;

    return  ( ( ($lambda ** $k) * exp(-$lambda) ) / fac($k) );
}

sub uniqueKmerThreshold($$$$) {
    my $wrk       = shift @_;
    my $asm       = shift @_;
    my $merSize   = shift @_;
    my $loss      = shift @_;
    my $bin       = getBinDirectory();
    my $errorRate = estimateRawError($wrk, $asm, "cor", $merSize);

    my $readLength = getNumberOfBasesInStore($wrk, $asm) / getNumberOfReadsInStore ($wrk, $asm);
    my $effective_coverage = getExpectedCoverage($wrk, $asm) * ( ($readLength - $merSize + 1)/$readLength ) * (1 - $errorRate) ** $merSize;

    my $threshold = 0;
    my $kMer_loss = poisson_pdf($effective_coverage, 0);

    return 1 if($kMer_loss > $loss); 

    my $keepTrying = 1;
    while($keepTrying)
    {
        $keepTrying = 0;
        my $p_true_kMers_threshold_p1 = poisson_pdf($effective_coverage, $threshold+1);
        if(($kMer_loss + $p_true_kMers_threshold_p1) <= $loss)
        {
                $threshold++;
                $kMer_loss += $p_true_kMers_threshold_p1;
                $keepTrying = 1;
        }
   }

   return ($threshold == 0 ? 1 : $threshold);
}

sub computeSampleSize($$$$$) {
    my $wrk      = shift @_;
    my $asm      = shift @_;
    my $tag      = shift @_;
    my $percent  = shift @_;
    my $coverage = shift @_;
    my $sampleSize = 0;

    my $minSampleSize = 100;
    my $maxSampleSize = getGlobal("${tag}MhapBlockSize") * 4;

   if (defined($percent)) {
      $sampleSize = int($percent * getNumberOfReadsInStore ($wrk, $asm))+1;
      $sampleSize++ if ($sampleSize % 2 != 0);
   } elsif (defined($coverage)) {
      $sampleSize = int(($coverage * getGlobal("genomeSize")) / (getNumberOfBasesInStore($wrk, $asm) / getNumberOfReadsInStore ($wrk, $asm))) + 1;
   }

   $sampleSize = $maxSampleSize if (defined($percent) && $sampleSize > $maxSampleSize);
   return $sampleSize < $minSampleSize ? $minSampleSize : $sampleSize;
}

sub runMHAP($$$$$$$$$$$$) {
    my ($wrk, $tag, $numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer, $sampleSize, $hash, $query, $out, $err) = @_;

    my $filterThreshold = getGlobal("${tag}MhapFilterThreshold");
    my $merSize         = getGlobal("${tag}MhapMerSize");
    my $javaPath        = getGlobal("java");
    my $bin             = getBinDirectory();

    print STDERR "--\n";
    print STDERR "-- PARAMETERS: hashes=$numHashes, minMatches=$minNumMatches, threshold=$threshold\n";
    print STDERR "--\n";   

    my $cmd  = "$javaPath -d64 -server -Xmx4g -jar $bin/mhap-" . getGlobal("${tag}MhapVersion") . ".jar ";
    $cmd .= "  --no-self --repeat-weight 0.9 -k $merSize --num-hashes $numHashes --num-min-matches $minNumMatches --ordered-sketch-size $ordSketch --ordered-kmer-size $ordSketchMer  --threshold $threshold --filter-threshold $filterThreshold --num-threads " . getGlobal("${tag}mhapThreads");
    $cmd .= " -s $hash -q $query  2> /dev/null | awk '{if (\$1 != \$2+$sampleSize) { print \$0}}' | $bin/errorEstimate -d 2 -m 0.95 -S - > $out 2> $err";
    runCommand($wrk, $cmd);
}


sub estimateRawError($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $merSize = shift @_;
    my $bin     = getBinDirectory();
    my $numReads = getNumberOfReadsInStore ($wrk, $asm);

    goto allDone   if (skipStage($WRK, $asm, "errorEstimate") == 1);
    goto allDone   if (-e "$wrk/asm.gkpStore/raw.estimate.out");
    goto allDone   if (getGlobal("errorrate") > 0);

    my ($numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer);

    $numHashes          =  10000;
    $minNumMatches      =    3;
    $threshold          =    0.65;
    $ordSketch          = 10000;
    $ordSketchMer       = getGlobal("${tag}MhapOrderedMerSize");

    # subsample raw reads
    my $sampleSize = computeSampleSize($wrk, $asm, $tag, 0.01, undef);
    $sampleSize /= 2;
    my $cmd = "$bin/gatekeeperDumpFASTQ -G $wrk/$asm.gkpStore -nolibname -fasta -r 1-$sampleSize -o - > $wrk/$asm.gkpStore/subset.fasta 2> /dev/null";
    runCommandSilently($wrk, $cmd, 1);
    my $min = $numReads - $sampleSize + 1;
    my $cmd = "$bin/gatekeeperDumpFASTQ -G $wrk/$asm.gkpStore -nolibname -fasta -r $min-$numReads -o - >> $wrk/$asm.gkpStore/subset.fasta 2> /dev/null";
    runCommandSilently($wrk, $cmd, 1);
    my $querySize = computeSampleSize($wrk, $asm, $tag, undef, 2);
    my $cmd = "$bin/gatekeeperDumpFASTQ -G $wrk/$asm.gkpStore -nolibname -fasta -r 1-$querySize -o - > $wrk/$asm.gkpStore/reads.fasta 2> /dev/null";
    runCommandSilently($wrk, $cmd, 1);

    print STDERR "--\n";
    print STDERR "-- ESTIMATOR (mhap) (raw) (hash sample size=". ($sampleSize*2) . ") (query sample size=$querySize)\n";
    runMHAP($wrk, $tag, $numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer, $sampleSize*2, "$wrk/$asm.gkpStore/subset.fasta", "$wrk/$asm.gkpStore/reads.fasta", "$wrk/$asm.gkpStore/raw.estimate.out", "$wrk/$asm.gkpStore/raw.estimate.err");
    unlink("$wrk/$asm.gkpStore/subset.fasta");
    unlink("$wrk/$asm.gkpStore/reads.fasta");

  allDone:
    return 0.15 if (! -e "$wrk/$asm.gkpStore/raw.estimate.out");

    my $errorRate = 0;
    open(L, "< $wrk/$asm.gkpStore/raw.estimate.out") or caExit("can't open '$wrk/$asm.gkpStore/raw.estimate.out' for reading: $!", undef);
    while (<L>) {
        $errorRate = sprintf "%.3f", ($_ / 2);
        $errorRate = 0.15 if ($errorRate <= 0.005);
    }
    close(L);

    return $errorRate;
}

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

    #  Mhap parameters - filterThreshold needs to be a string, else it is printed as 5e-06.
    #

    my ($numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer);

    $numHashes     	=  256;
    $minNumMatches 	=    4;
    $threshold     	=    0.85;
    $ordSketch     	= 1000;
    $ordSketchMer  	= getGlobal("${tag}MhapOrderedMerSize") + 2;

    make_path("$path");

    # subsample corrected reads, this assumes the fasta records are on a single line. We take some reads from the top and bottom of file to avoid sampling one library
    my $sampleSize = computeSampleSize($wrk, $asm, $tag, 0.01, undef);
    my $cmd = "gunzip -c $WRK/asm.correctedReads.fasta.gz |head -n $sampleSize > $path/subset.fasta";
    runCommandSilently($path, $cmd, 1);
    my $cmd = "gunzip -c $WRK/asm.correctedReads.fasta.gz |tail -n $sampleSize >> $path/subset.fasta";
    runCommandSilently($path, $cmd, 1);
    my $querySize =  computeSampleSize($wrk, $asm, $tag, undef, 2);
    my $cmd = "gunzip -c $WRK/asm.correctedReads.fasta.gz |head -n $querySize > $path/reads.fasta";
    runCommandSilently($path, $cmd, 1);
    my $cmd = "gunzip -c $WRK/asm.correctedReads.fasta.gz |tail -n $querySize >> $path/reads.fasta";
    runCommandSilently($path, $cmd, 1);

    # now compute the overlaps
    print STDERR "--\n";
    print STDERR "-- ESTIMATOR (mhap) (corrected) (hash sample size=$sampleSize) (query sample size=$querySize)\n";

    runMHAP($wrk, $tag, $numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer, $sampleSize, "$path/subset.fasta", "$path/reads.fasta", "$path/$asm.estimate.out", "$path/$asm.estimate.err");
    unlink("$path/subset.fasta");
    unlink("$path/reads.fasta");
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
    } elsif ($errorRate < 0.005) {
        print STDERR "-- Estimated error rate: " . ($errorRate*100) . "%, increasing to " . (0.005 * 100). "%.\n";
        $errorRate = 0.005;
    } else {
       print STDERR "-- Estimated error rate: " . ($errorRate * 100) . "%.\n";
    }
    setErrorRate($errorRate);
    showErrorRates("--  ");
    print STDERR "-- \n";
}
