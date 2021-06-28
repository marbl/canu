
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

package canu::Unitig;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(reportUnitigSizes unitig unitigCheck);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);
use POSIX qw(ceil);

use canu::Defaults;
use canu::Execution;

use canu::Configure;    #  For displayGenomeSize
use canu::SequenceStore;
use canu::Meryl;
use canu::Report;

use canu::Grid_Cloud;



sub reportBestError ($) {
    my $asm       = shift @_;
    my $log       = "unitigging/4-unitigger/$asm.best.report";
    my $report    = "";

    fetchFile($log);

    if (-e "$log") {
        open(F, "< $log") or caExit("failed to open filterOverlaps.log for reading: $!", undef);
        while (<F>) {
            $report .= "--  $_";
        }
        close(F);
    }

    addToReport("error rates", $report);
}



sub reportUnitigSizes ($$$) {
    my $asm       = shift @_;
    my $version   = shift @_;
    my $label     = shift @_;

    my $bin       = getBinDirectory();
    my $cmd       = "";

    my $ctgNum    = 0;
    my $ctgBases  = 0;
    my $ctgSizes  = "";

    my $usmNum   = 0;
    my $usmBases = 0;

    my $bubNum   = 0;
    my $bubBases = 0;

    my $rptNum    = 0;
    my $rptBases  = 0;

    my $gs        = getGlobal("genomeSize");

    my $V = substr("000000" . $version, -3);
    my $N = "$asm.ctgStore/seqDB.v$V.sizes.txt";

    fetchFile("unitigging/$N");

    if (! -e "unitigging/$N") {
        fetchFile("unitigging/$asm.ctgStore/seqDB.v$V.dat");
        fetchFile("unitigging/$asm.ctgStore/seqDB.v$V.tig");

        if (-e "unitigging/$asm.ctgStore/seqDB.v$V.tig") {
            $cmd  = "$bin/tgStoreDump \\\n";                     #  Duplicated at the end of unitigger.sh
            $cmd .= "  -S ../$asm.seqStore \\\n";
            $cmd .= "  -T ./$asm.ctgStore $version \\\n";
            $cmd .= "  -sizes -s " . getGlobal("genomeSize") . " \\\n";
            $cmd .= "> ./$N";

            if (runCommand("unitigging", $cmd)) {
                caExit("failed to generate unitig sizes", undef);
            }

            stashFile("unitigging/$N");
        }
    }

    #  If bogart is told to stop early, there won't be a '$N', but it isn't
    #  an error.  Just skip the stats.

    if (-e "unitigging/$N") {
        $ctgSizes .= "--            NG (bp)  LG (contigs)    sum (bp)\n";
        $ctgSizes .= "--         ----------  ------------  ----------\n";

        open(F, "< unitigging/$N") or caExit("failed to open 'unitigging/$N' for reading: $!", undef);
        while (<F>) {
            $rptBases  = $1  if (m/lenSuggestRepeat\s+sum\s+(\d+)/);
            $rptNum    = $1  if (m/lenSuggestRepeat\s+num\s+(\d+)/);

            $usmBases  = $1  if (m/lenUnassembled\s+sum\s+(\d+)/);
            $usmNum    = $1  if (m/lenUnassembled\s+num\s+(\d+)/);

            $bubBases  = $1  if (m/lenBubble\s+sum\s+(\d+)/);
            $bubNum    = $1  if (m/lenBubble\s+num\s+(\d+)/);

            $ctgBases  = $1  if (m/lenContig\s+sum\s+(\d+)/);
            $ctgNum    = $1  if (m/lenContig\s+num\s+(\d+)/);

            if (m/lenContig\s+ng(\d+)\s+(\d+)\s+bp\s+lg(\d+)\s+(\d+)\s+sum\s+(\d+)\s+bp$/) {
                my $n  = substr("          $1", -3);
                my $ng = substr("          $2", -10);
                my $lg = substr("          $4", -6);
                my $ss = substr("          $5", -10);

                $ctgSizes .= "--    $n  $ng        $lg  $ss\n";
            }

            #  User could have changed genome size since this report was dumped.
            if (m/lenContig.*genomeSize\s+(\d+)/) {
                $gs = $1;
            }
        }
        close(F);

        my $report;

        $report .= "-- Found, in version $version, $label:\n";
        $report .= "--   contigs:      $ctgNum sequences, total length $ctgBases bp (including $rptNum repeats of total length $rptBases bp).\n";
        $report .= "--   bubbles:      $bubNum sequences, total length $bubBases bp.\n";
        $report .= "--   unassembled:  $usmNum sequences, total length $usmBases bp.\n";
        $report .= "--\n";
        $report .= "-- Contig sizes based on genome size " . displayGenomeSize($gs) . "bp:\n";
        $report .= "--\n";
        $report .= $ctgSizes;
        $report .= "--\n";

        #  Hmmm.  Report wants a small tag, but $label is a bit verbose.

        addToReport("contigs",   $report)    if ($label eq "after unitig construction");
        addToReport("consensus", $report)    if ($label eq "after consensus generation");
    }
}



sub unitig ($) {
    my $asm     = shift @_;
    my $path    = "unitigging/4-unitigger";

    goto allDone    if (fileExists("unitigging/$asm.ctgStore/seqDB.v001.tig"));

    make_path($path)  if (! -d $path);

    #  If this is the first time here, reset the iteration count.  If the
    #  script exists already, we're doing a retry and so don't want to reset
    #  the iteration.
    #
    #  Parameter changes will still work, since Canu must be restarted
    #  manually to change the parameters, and that resets the iteration count
    #  implicitly.

    if (!fileExists("$path/unitigger.sh")) {
        resetIteration("unitig");
    }

    #  Dump a script to run the unitigger.

    open(F, "> $path/unitigger.sh") or caExit("can't open '$path/unitigger.sh' for writing: $!\n", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F "\n";
    print F fetchSeqStoreShellCode($asm, $path, "");
    print F "\n";
    print F fetchOvlStoreShellCode($asm, $path, "");
    print F "\n";
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    print F "#  Check if the outputs exist.  If they do, just quit.  (The boilerplate\n";
    print F "#  function for doing this fails if the file isn't strictly below the\n";
    print F "#  current directory, so some gymnastics is needed.)\n";
    print F "\n";
    print F "cd ..\n";
    print F fileExistsShellCode("exists", "unitigging", "$asm.ctgStore/seqDB.v001.tig", "");
    print F "cd 4-unitigger\n";
    print F "\n";
    print F "if [ \$exists = true ] ; then\n";
    print F "  echo Unitigger outputs exist, stopping.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";

    print F "\$bin/bogart \\\n";
    print F "  -S ../../$asm.seqStore \\\n";
    print F "  -O    ../$asm.ovlStore \\\n";
    print F "  -o     ./$asm \\\n";
    print F "  -gs "          . getGlobal("genomeSize")         . " \\\n";
    print F "  -eg "          . getGlobal("utgErrorRate")       . " \\\n";
    print F "  -eM "          . getGlobal("utgErrorRate")       . " \\\n";
    print F "  -mo "          . getGlobal("minOverlapLength")   . " \\\n";
    print F "  -covgapolap "  . getGlobal("minOverlapLength")   . " \\\n";
    print F "  -covgaptype "  . getGlobal("utgChimeraType")     . " \\\n";
    print F "  -lopsided 25 "                                   . " \\\n";
    print F "  -minolappercent   0.0 "                          . " \\\n";
    print F "  -dg "          . getGlobal("utgGraphDeviation")  . " \\\n";
    print F "  -db "          . getGlobal("utgBubbleDeviation") . " \\\n";
    print F "  -dr "          . getGlobal("utgRepeatDeviation") . " \\\n";
    print F "  -ca "          . getGlobal("utgRepeatConfusedBP"). " \\\n";
    print F "  -cp "          . getGlobal("utgRepeatConfusedPC"). " \\\n";
    print F "  -threads "     . getGlobal("batThreads")         . " \\\n"   if (defined(getGlobal("batThreads")));
    print F "  -M "           . getGlobal("batMemory")          . " \\\n"   if (defined(getGlobal("batMemory")));
    print F "  -unassembled " . getGlobal("contigFilter")       . " \\\n"   if (defined(getGlobal("contigFilter")));
    print F "  "              . getGlobal("batOptions")         . " \\\n"   if (defined(getGlobal("batOptions")));
    print F "> ./unitigger.err 2>&1\n";

    print F "\n";
    print F "if [ $? -ne 0 ] ; then\n";
    print F "  echo Unitigger failed to complete successfully.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "#  Note that we have finished.  All the rest is optional; if\n";
    print F "#  bogart is told to stop early, there won't be a ctgStore either.\n";
    print F "\n";
    print F "touch ./unitigger.finished\n";
    print F "\n";
    print F "if [ -e ./$asm.ctgStore ] ; then\n";
    print F "  mv ./$asm.ctgStore ../$asm.ctgStore\n";
    print F "fi\n";
    print F "\n";
    print F "#\n";
    print F "#  Generate tig sizes here, so we can avoid another pull of\n";
    print F "#  seqStore and ctgStore on cloud runs.  The output is generated\n";
    print F "#  again if it is missing, so no damage done if it isn't made here.\n";
    print F "#\n";
    print F "if [ -e ../$asm.ctgStore -a \\\n";
    print F "   ! -e ../$asm.ctgStore/seqDB.v001.sizes.txt ] ; then\n";
    print F "  \$bin/tgStoreDump \\\n";                     #  Duplicated in reportUnitigSizes()
    print F "    -S ../../$asm.seqStore \\\n";              #  Done here so we don't need another
    print F "    -T ../$asm.ctgStore 1 \\\n";               #  pull of seqStore and ctgStore
    print F "    -sizes -s " . getGlobal("genomeSize") . " \\\n";
    print F "   > ../$asm.ctgStore/seqDB.v001.sizes.txt\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e ../$asm.ctgStore ] ; then\n";
    print F "  cd ../$asm.ctgStore\n";
    print F stashFileShellCode("unitigging/$asm.ctgStore", "seqDB.v001.dat", "  ");
    print F stashFileShellCode("unitigging/$asm.ctgStore", "seqDB.v001.tig", "  ");
    print F stashFileShellCode("unitigging/$asm.ctgStore", "seqDB.v001.sizes.txt", "  ");
    print F "  cd -\n";
    print F "fi\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

    makeExecutable("$path/unitigger.sh");
    stashFile("$path/unitigger.sh");

  finishStage:
    generateReport($asm);
    #resetIteration("unitig");   #  Do not reset iteration here as this function is called every loop.

  allDone:
}




sub unitigCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "unitigging/4-unitigger";

    goto allDone      if (fileExists("$path/unitigger.success"));
    goto finishStage  if (fileExists("$path/unitigger.finished"));
    goto finishStage  if (fileExists("unitigging/$asm.ctgStore/seqDB.v001.tig"));

    fetchFile("$path/unitigger.sh");

    #  Since there is only one job, if we get here, we're not done.  Any other 'check' function
    #  shows how to process multiple jobs.  This only checks for the existence of the final outputs.
    #  (meryl-process, unitig, overlapStoreSequential are the same)

    #  If too many attempts, give up.

    if ($attempt >= getGlobal("canuIterationMax")) {
        print STDERR "--\n";
        print STDERR "-- Bogart failed, tried $attempt times, giving up.\n";
        print STDERR "--\n";
        caExit(undef, "$path/unitigger.err")   if (-e "$path/unitigger.err");
        caExit(undef, undef);
    }

    if ($attempt > 0) {
        print STDERR "--\n";
        print STDERR "-- Bogart failed, retry\n";
        print STDERR "--\n";
    }

    #  Otherwise, run some jobs.

    generateReport($asm);

    submitOrRunParallelJob($asm, "bat", $path, "unitigger", (1));
    return;

  finishStage:
    print STDERR "-- Unitigger finished successfully.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/unitigger.success") or caExit("can't open '$path/unitigger.success' for writing: $!", undef);
    close(F);

    stashFile("$path/unitigger.success");

    reportBestError($asm);
    reportUnitigSizes($asm, 1, "after unitig construction");

    generateReport($asm);
    resetIteration("unitigCheck");

  allDone:
    stopAfter("unitig");
}
