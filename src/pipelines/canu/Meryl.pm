package canu::Meryl;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(meryl);

use strict;

use canu::Defaults;
use canu::Execution;


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


sub meryl ($$$) {
    my $WRK = shift @_;
    my $wrk = $WRK;
    my $asm = shift @_;
    my $tag = shift @_;
    my $bin = getBinDirectory();
    my $cmd;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");

    my ($merSize, $merThresh, $merScale, $merDistinct, $merTotal);
    my ($ffile, $ofile);

    if (getGlobal("${tag}Overlapper") eq "ovl") {
        $merSize      = getGlobal("${tag}OvlMerSize");
        $merThresh    = getGlobal("${tag}OvlMerThreshold");
        $merScale     = 1.0;
        $merDistinct  = getGlobal("${tag}OvlMerDistinct");
        $merTotal     = getGlobal("${tag}OvlMerTotal");

        $ffile = "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta";   #  The fasta file we should be creating.
        $ofile = "$wrk/0-mercounts/$asm.ms$merSize";                      #  The meryl database 'intermediate file'.

    } elsif (getGlobal("${tag}Overlapper") eq "mhap") {
        $merSize      = getGlobal("${tag}mhapMerSize");
        $merThresh    = undef;
        $merScale     = 1.0;
        $merDistinct  = undef;
        $merTotal     = undef;

        $ffile = "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore";  #  The mhap-specific file we should be creating.
        $ofile = "$wrk/0-mercounts/$asm.ms$merSize";                      #  The meryl database 'intermediate file'.

    } else {
        caFailure("unknown ${tag}Overlapper '" . getGlobal("${tag}Overlapper") . "'", undef);
    }

    #  If the frequent mer file exists, don't bother running meryl.  We don't really need the
    #  databases.

    my $numTotal    = 0;
    my $numDistinct = 0;
    my $numUnique   = 0;
    my $largest     = 0;

    if (-e "$ofile.histogram.info") {
        open(F, "< $ofile.histogram.info") or caFailure("can't open meryl histogram information file '$ofile.histogram.info' for reading: $!\n", undef);
        while (<F>) {
            $numTotal    = $1   if (m/Found\s(\d+)\s+mers./);
            $numDistinct = $1   if (m/Found\s(\d+)\s+distinct\smers./);
            $numUnique   = $1   if (m/Found\s(\d+)\s+unique\smers./);
            $largest     = $1   if (m/Largest\smercount\sis\s(\d+)/);
        }
        close(F);
    }

    if ($numTotal > 0) {
        print STDERR "--  Found $numTotal $merSize-mers; $numDistinct distinct and $numUnique unique.  Largest count $largest.\n";
    }

    goto stopAfter  if (skipStage($WRK, $asm, "$tag-meryl") == 1);        #  Finished.
    goto allDone    if (-e "$ffile");

    #  Make a work space.

    system("mkdir $wrk/0-mercounts")  if (! -d "$wrk/0-mercounts");

    #  User supplied mers?  Just symlink to them.

    if (defined(getGlobal("${tag}OvlFrequentMers"))) {
        my $ffile = "$wrk/0-mercounts/$asm.frequentMers.fasta";
        my $sfile = getGlobal("${tag}OvlFrequentMers");

        if (! -e $ffile) {
            caFailure("${tag}OvlFrequentMers '$sfile' not found", undef)  if (! -e $sfile);
            print STDERR "Using frequent mers in '$sfile'\n";
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
        print STDERR "Threshold zero.  Empty file.\n";
        touch($ffile);

        goto allDone;
    }

    #  Build the database.

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
                print STDERR "Supplied merThreshold $merThresh is the smallest.\n";
                last;
            }

            if ((defined($merDistinct)) && ($merDistinct <= $distinct)) {
                $merThresh = (($merThresh > 0) && ($merThresh < $threshold)) ? $merThresh : $threshold;
                print STDERR "Supplied merDistinct $merDistinct with threshold $threshold is the smallest.\n";
                last;
            }

            if ((defined($merTotal)) && ($merTotal <= $total)) {
                $merThresh = (($merThresh > 0) && ($merThresh < $threshold)) ? $merThresh : $threshold;
                print STDERR "Supplied merTotal $merTotal with threshold $threshold is the smallest.\n";
                last;
            }
        }
        close(F);
    }

    #  Plot the histogram - annotated with the thesholds

    if (! -e "$ofile.histogram.png") {
        open(F, "> $ofile.histogram.gp");
        print F "\n";
        print F "\n";
        print F "\n";
        print F "unset multiplot\n";
        print F "\n";
        print F "set terminal png size 1000,800\n";
        print F "set output \"$ofile.histogram.png\"\n";
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

        runCommandSilently("$wrk/0-mercounts", "gnuplot $ofile.histogram.gp > /dev/null 2>&1");
    }

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
        open(O, "> $ofile.frequentMers.ignore")              or die "Failed to open '$ofile.frequentMers.mhap_ignore' for writing: $!\n";

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
    }

    #  Report the new threshold.

    if ((getGlobal("${tag}Overlapper") eq "ovl") && ($merThresh > 0) && (getGlobal("${tag}OvlMerThreshold") ne $merThresh)) {
        print STDERR "Reset ${tag}OvlMerThreshold from ", getGlobal("${tag}OvlMerThreshold"), " to $merThresh.\n";
        setGlobal("${tag}OvlMerThreshold", $merThresh);
    }

  purgeIntermediates:
    unlink "$ofile.mcidx"  if (getGlobal("saveMerCounts") == 0);
    unlink "$ofile.mcdat"  if (getGlobal("saveMerCounts") == 0);

  allDone:
    emitStage($WRK, $asm, "$tag-meryl");
  stopAfter:
    stopAfter("meryl");
}
