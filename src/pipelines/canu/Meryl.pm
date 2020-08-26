
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

package canu::Meryl;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(merylConfigure merylCountCheck merylProcessCheck);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);
use File::Basename;
use File::Copy;
use POSIX qw(floor ceil);

use canu::Defaults;
use canu::Execution;

use canu::SequenceStore;
use canu::Report;

use canu::Grid_Cloud;



sub merylParameters ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;

    my ($base, $path, $name, $merSize);

    #  Find a place to run stuff.

    $base = "haplotype"   if ($tag eq "hap");
    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/0-mercounts";
    $path = "$path-$asm"  if ($tag eq "hap");

    #  Decide on which set of parameters we need to be using, and make output file names.

    if (getGlobal("${tag}Overlapper") eq "ovl") {
        $merSize = getGlobal("${tag}OvlMerSize");
        $name    = "$asm.ms$merSize";

    } elsif (getGlobal("${tag}Overlapper") eq "mhap") {
        $merSize = getGlobal("${tag}mhapMerSize");
        $name    = "$asm.ms$merSize";

    } elsif (getGlobal("${tag}Overlapper") eq "minimap") {
        $merSize = 0;
        $name    = undef;

    } else {
        caFailure("unknown ${tag}Overlapper '" . getGlobal("${tag}Overlapper") . "'", undef);
    }

    #  Return all this goodness.  Well, there used to be a whole lot more stuff here, but
    #  it was never used and/or obsoleted by meryl now dumping mhap more or less directly.

    return($base, $path, $name, $merSize);
}



sub merylGenerateHistogram ($$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $hist;

    my ($base, $path, $name, $merSize) = merylParameters($asm, $tag);

    fetchFile("$path/$name.histogram");

    return(undef)   if (! -e "$path/$name.histogram");

    #  Load histogram data, limited to a max count of 30,000, to prevent
    #  these arrays from exhausting memory on very very frequent kmers.  It's
    #  just for display anyway.

    my $numTotal    = 0;
    my $numDistinct = 0;
    my $numUnique   = 0;
    my $largest     = 0;

    my $maxCount;
    my @numDistinct;
    my @fractDistinct;
    my @fractTotal;

    open(F, "< $path/$name.histogram");
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        $numUnique   = $1   if (m/unique\s+(\d+)\s/);
        $numDistinct = $1   if (m/distinct\s+(\d+)\s/);
        $numTotal    = $1   if (m/present\s+(\d+)\s/);

        if (($v[0] =~ m/^\d+$/) && ($v[0] < 30000)) {
            $maxCount             = $v[0];  #  maxCount is last count seen; histogram is sorted.
            $numDistinct[$v[0]]   = $v[1];
            $fractDistinct[$v[0]] = $v[2];
            $fractTotal[$v[0]]    = $v[3];
        }
    }
    close(F);

    if ($numTotal == 0) {
        caExit("meryl histogram '$path/$name.histogram' shows no kmers", undef);
    }

    #  Find blocks for the histogram.

    my @TD;  #  Total number of distinct kmers in the i'th histogram block.
    my @FU;  #  Fraction distinct ...
    my @FT;  #  Fraction total    ...

    my $TDmax  = 0;   #  Max count of of any block, excluding the first (we ignore the tail of this block when drawing the histogram)

    my $maxRows = 40;   #  Number of rows for the count histogram.  We'll change the bin
    my $rowIncr = 1;    #  size (rowIncr) based on the depth of coverage.

    #  Find block sizes for the histogram.

    my $lo = 1;
    my $hi = 2;
    my $st = 1;

    for (my $bb=1; $bb<6; $bb++) {
        my $ft;

        $lo = 1;
        $hi = 2;
        $st = 1;
        for (my $ii=0; $ii <= $maxRows; $ii++) {
            for (my $jj=$lo; $jj < $hi; $jj++) {
                $ft = $fractTotal[$jj]   if (defined($fractTotal[$jj]))
            }

            $lo  = $hi;
            $hi += $st;
            $st += $bb;
        }

        if ($ft >= 0.95) {
            $rowIncr = $bb;
            last;
        }
    }

    #  Find blocks for the histogram.

    $lo = 1;
    $hi = 2;
    $st = 1;
    for (my $ii=0; $ii <= $maxRows; $ii++) {
        for (my $jj=$lo; $jj < $hi; $jj++) {
            $TD[$ii] += $numDistinct[$jj];                                                   #  Sum the number of distinct mers we've seen

            $FU[$ii] = ($fractDistinct[$ii] < $FU[$ii]) ? $FU[$ii] : $fractDistinct[$jj];    #  But the fractions are already cumulative,
            $FT[$ii] = ($fractTotal[$ii]    < $FT[$ii]) ? $FT[$ii] : $fractTotal[$jj];       #  we just need to skip zeros.
        }

        if ($ii > 0) {
            $TDmax = ($TDmax < $TD[$ii]) ? $TD[$ii] : $TDmax;
        }

        $lo  = $hi;
        $hi += $st;
        $st += $rowIncr;
    }

    if ($TDmax == 0) {           #  A pathological case; if all kmers are unique,
        $TDmax = $TD[0];         #  no max size is set, Xscale is zero,
    }                            #  and we fail.

    my $maxY   = 1;              #  Last count to include in the histogram.
    my $Xscale = $TDmax / 70;    #  Scale of each * in the histogram.

    for (my $ii=0; $ii <= $maxRows; $ii++) {
        $maxY = $ii  if ($TD[$ii] > 0);
    }

    #  Now just draw the histogram

    $hist .= "--\n";
    $hist .= "--  $merSize-mers                                                                                           Fraction\n";
    $hist .= "--    Occurrences   NumMers                                                                         Unique Total\n";

    $lo = 1;
    $hi = 2;
    $st = 1;
    for (my $ii=0; $ii<=$maxRows; $ii++) {
        my $numXs = int($TD[$ii] / $Xscale);

        if ($numXs <= 70) {
            $hist .= sprintf("--  %6d-%6d %9d %s%s %.4f %.4f\n",
                             $lo, $hi-1, $TD[$ii],
                             "*" x      ($numXs),
                             " " x (70 - $numXs), $FU[$ii], $FT[$ii]);
        } else {
            $hist .= sprintf("--  %6d-%6d %9d %s%s %.4f %.4f\n",
                             $lo, $hi-1, $TD[$ii],
                             "*" x 67,
                             "-->", $FU[$ii], $FT[$ii]);
        }

        last   if ($ii >= $maxY);

        $lo  = $hi;
        $hi += $st;
        $st += $rowIncr;
    }

    $hist .= sprintf("--\n");
    $hist .= sprintf("-- %11d (max occurrences)\n",             $largest);
    $hist .= sprintf("-- %11d (total mers, non-unique)\n",      $numTotal    - $numUnique);
    $hist .= sprintf("-- %11d (distinct mers, non-unique)\n",   $numDistinct - $numUnique);
    $hist .= sprintf("-- %11d (unique mers)\n",                 $numUnique);

    return($hist);
}



sub merylPlotHistogram ($$$$) {
    my $path   = shift @_;
    my $asm    = shift @_;
    my $name   = shift @_;
    my $size   = shift @_;  #  Size of image, not merSize!

    my $gnuplot = getGlobal("gnuplot");
    my $format  = getGlobal("gnuplotImageFormat");

    return  if (  fileExists("$path/$name.histogram.gp"));
    return  if (! defined($gnuplot) || ($gnuplot eq ""));
    return  if (! fileExists("$path/$name.histogram"));

    fetchFile("$path/$name.histogram");

    open(F, "> $path/$name.histogram.gp");
    print F "\n";
    print F "unset multiplot\n";
    print F "\n";
    print F "set terminal $format size $size,$size\n";
    print F "set output '$name.histogram.$format'\n";
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
    print F "plot [0.5:1.0] '$name.histogram' using 3:4 with lines title 'Distinct-vs-Total'\n";
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
    print F "plot [0.975:1.0] '$name.histogram' using 3:4 with lines title 'Distinct-vs-Total'\n";
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
    print F "plot [0:200] '$name.histogram' using 1:2 with lines title 'Histogram'\n";
    close(F);

    if (runCommandSilently($path, "$gnuplot < /dev/null ./$name.histogram.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

    stashFile("$path/$name.histogram.gp");
    stashFile("$path/$name.histogram.$format");
}



sub merylConfigure ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    my ($base, $path, $name, $merSize) = merylParameters($asm, $tag);

    goto allDone   if ($merSize == 0);

    goto allDone   if (fileExists("$path/meryl-count.sh") &&            #  Configure 'output' exists.
                       fileExists("$path/meryl-process.sh"));           #
    goto allDone   if (fileExists("$path/$name.dump") &&                #  Meryl output(s) exist
                       fileExists("$path/$name.ignore.gz"));            #

    make_path($path)  if (! -d $path);

    #
    #  User supplied mers?  Copy them to the proper location and exit.
    #

    if (defined(getGlobal("${tag}OvlFrequentMers"))) {
        my $merFile = getGlobal("${tag}OvlFrequentMers");

        caFailure("${tag}OvlFrequentMers '$merFile' not found", undef)  if (! -e $merFile);

        if (${tag} eq "cor") {
            touch("$path/$name.dump");
            copy($merFile, "$path/$name.ignore.gz");
        }

        else {
            copy($merFile, "$path/$name.dump");
            touch("$path/$name.ignore.gz");
        }

        stashFile("$path/$name.dump");
        stashFile("$path/$name.ignore.gz");

        goto allDone;
    }

    #
    #  Based on the number of reads in the store, decide on a preliminary set
    #  of batch sizes to use.
    #

    my $nr  = getNumberOfReadsInStore($asm, $tag);
    my $nb  = getNumberOfBasesInStore($asm, $tag);

    my $maxSplit = int(floor($nb / 100000000));

    $maxSplit = 1    if ($maxSplit == 0);
    $maxSplit = $nr  if ($nr < $maxSplit);

    #
    #  Let meryl run to decide how much memory it wants to use.
    #
    #  This is a pretty stupid way to do it; this really should be built into meryl itself.
    #

    my $mem = getGlobal("merylMemory");
    my $thr = getGlobal("merylThreads");
    my $cov = getExpectedCoverage($asm, $tag);

    open(F, "> $path/meryl-configure.sh");
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";                                     #  setWDSC() should be OK here, but issue
    #print F setWorkDirectoryShellCode($path);        #  #1740 had a problem with it, so it is
    print F fetchSeqStoreShellCode($asm, $path, "");  #  disabled.  We're not sure why it broke.

    foreach my $ss (qw(01 02 04 06 08 12 16 20 24 32 40 48 56 64 96)) {
        next  if ($ss > $maxSplit);

        print F "\n";
        print F "$bin/meryl -C k=$merSize threads=$thr memory=$mem \\\n";
        print F "  count segment=1/$ss ../../$asm.seqStore \\\n";
        print F "> $name.config.$ss.out 2>&1";
    }
    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/meryl-configure.sh");
    stashFile("$path/meryl-configure.sh");

    my $runConfigure = 1;

    foreach my $ss (qw(01 02 04 06 08 12 16 20 24 32 40 48 56 64 96)) {
        if (-e "$path/$name.config.$ss.out") {
            $runConfigure = 0;
            next;
        }
    }

    if ($runConfigure) {
        if (runCommand($path, "./meryl-configure.sh > ./meryl-configure.err 2>&1")) {
            caFailure("meryl failed to configure", "$path/meryl-configure.err");
        }
        unlink("$path/meryl-configure.err");
    }

    #
    #  Make sense of all the configuration attempts.  Pick the smallest number of segments that
    #  results in the smallest number of batches.
    #

    my $merylMemory;
    my $merylSegments;
    my $merylBatches  = 1048576;

    printf(STDERR "--  segments   memory batches\n");
    printf(STDERR "--  -------- -------- -------\n");

    foreach my $ss (qw(01 02 04 06 08 12 16 20 24 32 40 48 56 64 96)) {
        next  if ($ss > $maxSplit);

        my $mem = undef;
        my $bat = undef;

        if (! -e "$path/$name.config.$ss.out") {
            next;
        }

        #  This message comes from meryl/merylOp-count.C reportNumberOfOutputs().

        open(F, "< $path/$name.config.$ss.out") or caExit("can't open '$path/$name.config.$ss.out' for reading: $!", undef);
        while (<F>) {
            if (m/Configured\s+\w+\s+mode\s+for\s+(\d*.\d*)\s+GB\s+memory\s+per\s+batch,\s+and\s+up\s+to\s+(\d+)\s+batch/) {
                $mem = $1;
                $bat = $2;
            }
        }
        close(F);

        caExit("failed to parse meryl configure output '$path/$name.config.$ss.out'", "$path/$name.config.$ss.out")   if (!defined($mem) || !defined($bat));

        if ($bat < $merylBatches) {
            $merylMemory   = int($mem) + 2;
            $merylSegments =     $ss;
            $merylBatches  =     $bat;
        }

        printf(STDERR "--       %3s %5.2f GB     %3d\n", $ss, $mem, $bat);
    }

    print STDERR "--\n";
    print STDERR "--  For $nr reads with $nb bases, limit to $maxSplit batch", ($maxSplit == 1) ? "" : "es", ".\n";
    print STDERR "--  Will count kmers using $merylSegments jobs, each using $merylMemory GB and $thr threads.\n";
    print STDERR "--\n";

    setGlobal("merylMemory",  $merylMemory);
    setGlobal("merylThreads", $thr);    #  Redundant.

    #
    #  Build a script for running meryl.  A similar version is used in HaplotypeReads.pm.
    #

    open(F, "> $path/meryl-count.sh") or caExit("can't open '$path/meryl-count.sh' for writing: $!", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    print F "if [ \$jobid -gt $merylSegments ]; then\n";
    print F "  echo Error: Only $merylSegments jobs, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %02d \$jobid`\n";
    print F "\n";
    print F "#  If the meryl database exists, we're done.\n";
    print F "\n";
    print F "if [ -e ./$asm.\$jobid.meryl/merylIndex ] ; then\n";
    print F "  echo Kmers for batch \$jobid exist.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "#  If the meryl output exists in the object store, we're also done.\n";
    print F "\n";
    print F fileExistsShellCode("exist1", "$path", "$asm.\$jobid.meryl.tar.gz");
    print F "if [ \$exist1 = true ] ; then\n";
    print F "  echo Kmers for batch \$jobid exist in the object store.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "#  Nope, not done.  Fetch the sequence store.\n";
    print F "\n";
    print F fetchSeqStoreShellCode($asm, $path, "");
    print F "\n";
    print F "#  And compute.\n";
    print F "\n";
    print F "$bin/meryl k=$merSize threads=$thr memory=$merylMemory \\\n";
    print F "  count \\\n";
    print F "    segment=\$jobid/$merylSegments ../../$asm.seqStore \\\n";
    print F "    output ./$asm.\$jobid.meryl.WORKING \\\n";
    print F "&& \\\n";
    print F "mv -f ./$asm.\$jobid.meryl.WORKING ./$asm.\$jobid.meryl\n";
    print F "\n";

    if (defined(getGlobal("objectStore"))) {
        print F stashMerylShellCode($path, "$asm.\$jobid.meryl", "");
    }

    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/meryl-count.sh");
    stashFile("$path/meryl-count.sh");

    #
    #  Build a script for processing meryl results
    #

    my @jobs;

    for (my $seg=1; $seg<=$merylSegments; $seg++) {
        push @jobs, substr("00$seg", -2);
    }

    open(F, "> $path/meryl-process.sh") or caExit("can't open '$path/meryl-process.sh' for writing: $!", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    print F "if [ \$jobid -gt 1 ]; then\n";
    print F "  echo Error: Only 1 job, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";

    print F "#  If the meryl ignore files exst, then we're done.\n";
    print F "\n";
    print F "if [ -e ./$name.histogram -a -e ./$name.dump -a -e ./$name.ignore.gz ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "#  If those exist in the object store, we're also done.\n";
    print F "\n";
    print F fileExistsShellCode("exists1", "$path", "$name.histogram");
    print F fileExistsShellCode("exists2", "$path", "$name.dump");
    print F fileExistsShellCode("exists3", "$path", "$name.ignore.gz");
    print F "if [ \$exists1 = true -a \$exists2 = true -a \$exists3 = true ] ; then\n";
    print F "  echo \"Output files '$name.histogram', '$name.dump' and '$name.ignore.gz' exist in '$path'.\"\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "\n";
    print F "#  Nope, not done.  Fetch all the intermediate meryl databases.\n";
    print F "\n";

    if (defined(getGlobal("objectStore"))) {
        print F fetchMerylShellCode($path, "$asm.$_.meryl", "")   foreach (@jobs);   #  One line, yay, but not use of $_.
    }

    print F "\n";
    print F "#\n";
    print F "#  Merge counting jobs, strip out unique kmers.\n";
    print F "#\n";
    print F "\n";
    print F "if [ ! -e ./$name/merylIndex ] ; then\n";
    print F "  $bin/meryl threads=$thr memory=$merylMemory \\\n";
    print F "    greater-than 1 \\\n";
    print F "      output $name.WORKING \\\n";
    print F "      union-sum  \\\n";
    print F "        ./$asm.$_.meryl \\\n"   foreach (@jobs);   #  One line, yay, but not use of $_.
    print F "  && \\\n";
    print F "  mv -f ./$name.WORKING ./$name\n";
    print F "\n";
    print F "  #  Fail if there is no meryl database.\n";
    print F "  if [ ! -e ./$name/merylIndex ] ; then\n";
    print F "    echo meryl merge failed.\n";
    print F "    exit 1\n";
    print F "  fi\n";
    print F "\n";
    print F "  #  Remove meryl intermediate files.\n";
    print F "  rm -rf ./$asm.$_.meryl ./$asm.$_.meryl.err\n"       foreach (@jobs);   #  One line, yay, but not use of $_.
    print F "fi\n";
    print F "\n";
    print F "#\n";
    print F "#  Dump a histogram, 'cause they're useful.\n";
    print F "#\n";
    print F "\n";
    print F "if [ ! -e ./$name.histogram ] ; then\n";
    print F "  $bin/meryl threads=1 memory=1 \\\n";
    print F "    statistics ./$name \\\n";
    print F "  > ./$name.histogram.WORKING \\\n";
    print F "  && \\\n";
    print F "  mv ./$name.histogram.WORKING ./$name.histogram\n";
    print F "fi\n";
    print F "\n";

    my $mthresh   = undef;
    my $mdistinct = undef;
    my $mwordfreq = undef;

    if (getGlobal("${tag}Overlapper") eq "ovl") {
        $mthresh   = getGlobal("${tag}OvlMerThreshold");       #  Kmer must meet at least BOTH thresholds.
        $mdistinct = getGlobal("${tag}OvlMerDistinct");
        $mwordfreq = undef;
    }

    if (getGlobal("${tag}Overlapper") eq "mhap") {
        $mthresh   = int(5.0 * getExpectedCoverage($asm, $tag));
        $mdistinct = undef;
        $mwordfreq = getGlobal("${tag}MhapFilterThreshold");
    }

    if ("${tag}Overlapper" eq "minimap") {
    }

    print F "#\n";
    print F "#  Dump frequent mers.\n";
    print F "#\n";
    print F "#  The indenting of the at-least options is misleading.  'print'\n";
    print F "#  takes input from the first 'at-least', which that takes input from\n";
    print F "#  the second 'at-least'.  The effect is the same as taking the\n";
    print F "#  'intersection' of all the 'at-least' filters -- logically, it is\n";
    print F "#  doing 'at-least X AND at-least Y AND at-least Z'.\n";
    print F "#\n";
    print F "\n";
    print F "if [ ! -e ./$name.dump ] ; then\n";
    print F "  $bin/meryl threads=$thr memory=$merylMemory \\\n";
    print F "    print ./$name.##.dump \\\n";
    print F "      at-least distinct=$mdistinct \\\n"         if (defined($mdistinct));
    print F "      at-least threshold=$mthresh \\\n"          if (defined($mthresh));
    print F "      at-least word-frequency=$mwordfreq \\\n"   if (defined($mwordfreq));
    print F "        ./$name\n";
    print F "\n";
    print F "  cat ./$name.??.dump > ./$name.dump\n";
    print F "  rm -f ./$name.??.dump\n";
    print F "fi\n";
    print F "\n";
    print F "#\n";
    print F "#  Convert the dumped kmers into a mhap ignore list.\n";
    print F "#\n";
    print F "#    numKmers - number of kmers we're filtering\n";
    print F "#    totKmers - total number of kmers in the dataset\n";
    print F "\n";
    print F "if [ ! -e ./$name.ignore.gz ] ; then\n";
    print F "  numKmers=`wc -l < ./$name.dump`\n";
    print F "  totKmers=`$bin/meryl statistics ./$name | grep present | awk '{ print \$2 }'`\n";
    print F "\n";

    if (defined(getGlobal("objectStore"))) {
        print F fetchFileShellCode($path, "meryl-make-ignore.pl", "  ");
    }

    print F "\n";
    print F "  ./meryl-make-ignore.pl \$numKmers \$totKmers < ./$name.dump | gzip -1c > ./$name.ignore.gz\n";
    print F "fi\n";
    print F "\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Save the final meryl database.\n";
        print F "\n";
        print F stashMerylShellCode($path, $name, "");

        print F "\n";
        print F "#  Save the histogram.\n";
        print F "\n";
        print F stashFileShellCode($path, "$name.histogram", "");

        print F "\n";
        print F "#  Save the overlapInCore ignore file.\n";
        print F "\n";
        print F stashFileShellCode($path, "$name.dump", "");

        print F "\n";
        print F "#  Save the mhap ignore file.\n";
        print F "\n";
        print F stashFileShellCode($path, "$name.ignore.gz", "");

    }

    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/meryl-process.sh");
    stashFile("$path/meryl-process.sh");

    #
    #  Build a (tiny) script for converting meryl dumps into mhap ignore files.
    #

    open(F, "> $path/meryl-make-ignore.pl") or caExit("can't open '$path/meryl-make-ignore.pl' for writing: $!", undef);
    print F "#!/usr/bin/env perl\n";
    print F "\n";
    print F "my \$numKmers = shift @ARGV;  #  The number of kmers in the input dump file.\n";
    print F "my \$totKmers = shift @ARGV;  #  The total number of kmers in the reads.\n";
    print F "\n";
    print F "printf \"%s\\t%d\\n\", 0, \$numKmers * 2;\n";
    print F "\n";
    print F "while (<STDIN>) {\n";
    print F "  chomp;\n";
    print F "\n";
    print F "  my \@v = split '\\s+', \$_;\n";
    print F "\n";
    print F "  printf \"%s\\t%.10f\\t%d\\t%d\\n\", \$v[0], \$v[1] / \$totKmers, \$v[1], \$totKmers;\n";
    print F "\n";
    print F "  \$v[0] =  reverse \$v[0];\n";
    print F "  \$v[0] =~ tr/ACGT/TGCA/;\n";
    print F "\n";
    print F "  printf \"%s\\t%.10f\\t%d\\t%d\\n\", \$v[0], \$v[1] / \$totKmers, \$v[1], \$totKmers;\n";
    print F "}\n";
    print F "\n";
    print F "exit(0);\n";
    close(F);

    makeExecutable("$path/meryl-make-ignore.pl");
    stashFile("$path/meryl-make-ignore.pl");


  finishStage:
    generateReport($asm);
    resetIteration("merylConfigure");

  allDone:
    stopAfter("meryl-configure");
}



sub merylCountCheck ($$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $bin     = getBinDirectory();
    my $cmd;

    my ($base, $path, $name, $merSize) = merylParameters($asm, $tag);

    #  If the frequent mer file exists, don't bother running meryl.  We don't really need the
    #  databases.

    goto allDone      if ($merSize == 0);

    goto allDone      if (fileExists("$path/meryl-count.success"));
    goto allDone      if (fileExists("$path/$name.dump") &&
                          fileExists("$path/$name.ignore.gz"));

    goto finishStage  if (fileExists("$path/$name/merylIndex"));

    fetchFile("$path/meryl-count.sh");

    #  Scan the script to determine how many jobs there are.

    my $jobs = 0;

    open(F, "< $path/meryl-count.sh") or caExit("can't open '$path/meryl-count.sh' for reading: $!", undef);
    while (<F>) {
        if (m/Only\s(\d+)\sjobs/) {
            $jobs = $1;
        }
    }
    close(F);

    caExit("failed to find the number of jobs in '$path/meryl-count.sh'", undef)  if ($jobs == 0);

    #  Figure out if all the tasks finished correctly.

    my $currentJobID = "01";

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if      ((fileExists("$path/$asm.$currentJobID.meryl")) ||
                 (fileExists("$path/$asm.$currentJobID.meryl.tar.gz"))) {
            push @successJobs, "$path/$asm.$currentJobID.meryl\n";

        } else {
            $failureMessage .= "--   job $asm.$currentJobID.meryl FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Kmer counting (meryl-count) jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Kmer counting (meryl-count) jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "meryl", $path, "meryl-count", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " Kmer counting (meryl) outputs.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/meryl-count.success") or caExit("can't open '$path/meryl-count.success' for writing: $!", undef);
    close(F);

    stashFile("$path/meryl-count.success");

    generateReport($asm);
    resetIteration("$tag-merylCountCheck");

  allDone:
    stopAfter("meryl-count");
}



sub merylProcessCheck ($$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $bin     = getBinDirectory();
    my $cmd;

    my ($base, $path, $name, $merSize) = merylParameters($asm, $tag);

    #  If the frequent mer file exists, don't bother running meryl.  We don't really need the
    #  databases.

    goto allDone      if ($merSize == 0);

    goto allDone      if (fileExists("$path/meryl-process.success"));

    goto finishStage  if (fileExists("$path/$name.histogram") &&
                          fileExists("$path/$name.dump") &&
                          fileExists("$path/$name.ignore.gz"));

    fetchFile("$path/meryl-process.sh");

    #  Since there is only one job, if we get here, we're not done.  Any other 'check' function
    #  shows how to process multiple jobs.  This only checks for the existence of the final outputs.
    #  (alignGFA and unitig are the same)

    #  If too many attempts, give up.

    if ($attempt >= getGlobal("canuIterationMax")) {
        print STDERR "--\n";
        print STDERR "-- meryl-process failed, tried $attempt times, giving up.\n";
        print STDERR "--\n";
        caExit(undef, undef);
    }

    if ($attempt > 0) {
        print STDERR "--\n";
        print STDERR "-- meryl-process failed, retry.\n";
        print STDERR "--\n";
    }

    #  Otherwise, run some jobs.

    generateReport($asm);

    submitOrRunParallelJob($asm, "meryl", $path, "meryl-process", (1));
    return;

  finishStage:
    print STDERR "-- Meryl finished successfully.  Kmer frequency histogram:\n";

    merylPlotHistogram($path, $asm, $name, 1024);

    addToReport("${tag}Meryl", merylGenerateHistogram($asm, $tag));

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/meryl-process.success") or caExit("can't open '$path/meryl-process.success' for writing: $!", undef);
    close(F);

    stashFile("$path/meryl-process.success");

    generateReport($asm);
    resetIteration("meryl-process");

    if (getGlobal("saveMerCounts") == 0) {
        print STDERR "--\n";
        print STDERR "-- Removing meryl database '$path/$name'.\n";
        remove_tree("$path/$name")
    } else {
        print STDERR "--\n";
        print STDERR "-- Meryl database '$path/$name' saved because 'saveMerCounts=true'.\n";
    }

  allDone:
    stopAfter("meryl-process");
    stopAfter("meryl");
}
