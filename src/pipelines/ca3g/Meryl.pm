package ca3g::Meryl;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(meryl);

use strict;

use ca3g::Defaults;
use ca3g::Execution;


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


#  Generates   $wrk/0-mercounts/$asm-ms$merSize.frequentMers.fasta
#  stopBefore  meryl (stops before meryl itself runs)
#  stopAfter   meryl (stops after output is generated, even if it is just a symlink)


sub meryl ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;

    my $merSize      = getGlobal('ovlMerSize');
    my $merThresh    = getGlobal("ovlMerThreshold");
    my $merScale     = 1.0;
    my $merDistinct  = getGlobal("ovlMerDistinct");
    my $merTotal     = getGlobal("ovlMerTotal");

    my $ffile = "$wrk/0-mercounts/$asm-ms$merSize.frequentMers.fasta";   #  The fasta file we should be creating.
    my $ofile = "$wrk/0-mercounts/$asm-ms$merSize";                      #  The merl database 'intermediate file'.

    #  If the frequent mer file exists, don't bother running meryl.  We don't really need the
    #  databases.

    return  if (-e "$ffile");

    #  Make a work space.

    system("mkdir $wrk/0-mercounts")  if (! -d "$wrk/0-mercounts");

    #  User supplied mers?  Just symlink to them.

    if (defined(getGlobal("ovlFrequentMers"))) {
        my $ffile = "$wrk/0-mercounts/$asm.frequentMers.fasta";
        my $sfile = getGlobal("ovlFrequentMers");

        if (! -e $ffile) {
            caFailure("frequentMers '$sfile' not found", undef)  if (! -e $sfile);
            print STDERR "Using frequent mers in '$sfile'\n";
            symlink $sfile, $ffile;
        }

        stopAfter("meryl");
        return;
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

    print "$merThresh $merDistinct $merTotal\n";

    if ((defined($merThresh))    &&
        ($merThresh ne "auto")   &&
        ($merThresh == 0)        &&
        (!defined($merDistinct)) &&
        (!defined($merTotal))) {
        print STDERR "Threshold zero.  Empty file.\n";
        touch($ffile);
        stopAfter("meryl");
        return;
    }

    #  Build the database.

    if (! -e "$ofile.mcdat") {
        my $merylMemory  = getGlobal("merylMemory");
        my $merylThreads = getGlobal("merylThreads");

        if ($merylMemory !~ m/^-/) {
            $merylMemory = "-memory $merylMemory";
        }

        $cmd  = "$bin/meryl ";
        $cmd .= " -B -C -v -m $merSize $merylMemory -threads $merylThreads ";
        $cmd .= " -L 2 ";
        $cmd .= " -s $wrk/$asm.gkpStore ";
        $cmd .= " -o $ofile ";
        $cmd .= "> $wrk/0-mercounts/meryl.err 2>&1";

        stopBefore("meryl", $cmd);

        if (runCommand("$wrk/0-mercounts", $cmd)) {
            caFailure("meryl failed", "$wrk/0-mercounts/meryl.err");
        }
        unlink "$wrk/0-mercounts/meryl.err";
    }

    #  Dump a histogram.

    if (! -e "$ofile.histogram") {
        $cmd  = "$bin/meryl -Dh -s $ofile > $ofile.histogram 2> $ofile.histogram.err";

        if (runCommand("$wrk/0-mercounts", $cmd)) {
            rename "$ofile.histogram", "$ofile.histogram.FAILED";
            caFailure("meryl histogram failed", "$ofile.histogram.err");
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

    print STDERR "MERYL SHOULD BE PLOTTING A HISTOGRAM.\n";

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

    runCommand("$wrk/0-mercounts", "gnuplot $ofile.histogram.gp");


    #  Generate the frequent mers

    if (! -e $ffile) {
        $cmd  = "$bin/meryl -Dt -n $merThresh -s $ofile > $ffile 2> $ffile.err";

        if (runCommand("$wrk/0-mercounts", $cmd)) {
            unlink $ffile;
            caFailure("meryl failed to dump frequent mers", "$ffile.err");
        }

        unlink "$ffile.err";
    }

    #  Report the new threshold.

    if (($merThresh > 0) && (getGlobal("ovlMerThreshold") ne $merThresh)) {
        print STDERR "Reset ovlMerThreshold from ", getGlobal("ovlMerThreshold"), " to $merThresh.\n";
        setGlobal("ovlMerThreshold", $merThresh);
    }

    stopAfter("meryl");
}
