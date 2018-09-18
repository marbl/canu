
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
 #    src/pipelines/canu/CorrectReads.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-APR-09 to 2015-SEP-03
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-19
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-17
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::HaplotypeReads;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(haplotypeReadsExist
             haplotypeCountConfigure
             haplotypeCountCheck
             haplotypeMergeCheck
             haplotypeSubtractCheck
             haplotypeReadsConfigure
             haplotypeReadsCheck);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;

use canu::Configure;
use canu::SequenceStore;
use canu::Report;

use canu::Grid_Cloud;

use POSIX qw(ceil);



sub haplotypeReadsExist ($@) {
    my $asm            = shift @_;
    my @haplotypes     =       @_;

    #  If no haplotypes, no haplotype reads exist.  Of note, this string causes the canu
    #  main to skip all haplotype processing.

    if (scalar(@haplotypes) == 0){
        return("no-haplotypes");
    }

    #  Check for each file, failing if any one doesn't exist.

    foreach my $haplotype (@haplotypes) {
        return("no")   if (! -e "./haplotype/haplotype-$haplotype.fasta.gz");
    }

    return("yes");
}



sub haplotypeCountConfigure ($%) {
    my $asm            = shift @_;
    my %haplotypeReads =       @_;
    my $bin            = getBinDirectory();
    my $cmd;
    my $path           = "haplotype/0-kmers";

    my @haplotypes     = keys %haplotypeReads;

    goto allDone   if (skipStage($asm, "haplotypeCountConfigure") == 1);
    goto allDone   if (fileExists("$path/meryl-count.sh") &&
                       fileExists("$path/meryl-merge.sh") &&
                       fileExists("$path/meryl-subtract.sh"));

    make_path($path)  if (! -d $path);

    #
    #  For each haplotype, split the sequence files into fixed-size chunks.
    #

    foreach my $haplotype (@haplotypes) {
        my @readFiles     = split '\0', $haplotypeReads{$haplotype};
        my $fileNumber    = "001";
        my $fileLength    = 0;
        my $fileLengthMax = 100000000;

        next  if (-e "$path/reads-$haplotype.success");

        print STDERR "--\n";
        print STDERR "--  Split haplotype reads for '$haplotype' into multiple files.\n";
        print STDERR "--   --> '$path/reads-$haplotype-$fileNumber.fasta.gz'.\n";
        open(OUT, "| gzip -1c > $path/reads-$haplotype-$fileNumber.fasta.gz");

        foreach my $file (@readFiles) {  
            print STDERR "--   <-- '$file'.\n";
            open(INP, "< $file");
            open(INP, "gzip  -dc $file |")   if ($file =~ m/\.gz$/);
            open(INP, "bzip2 -dc $file |")   if ($file =~ m/\.ba2$/);
            open(INP, "xz    -dc $file |")   if ($file =~ m/\.xz$/);

            while (!eof(INP)) {
                $_ = <INP>;

                #  If looks like FASTA or FASTQ name, make new sequence
                if ((m/^\@/) ||
                    (m/^>/)) {
                    print OUT ">\n";
                    next;
                }

                #  If looks like FASTQ quality values, skip them
                if (m/^\+/) {
                    $_ = <INP>;
                    next;
                }

                #  Otherwise, assume it's sequence to keep, print it
                print OUT $_;

                $fileLength += length($_) - 1;

                if ($fileLength > $fileLengthMax) {
                    close(OUT);

                    $fileNumber++;
                    $fileLength = 0;

                    print STDERR "--   --> '$path/reads-$haplotype-$fileNumber.fasta.gz'.\n";
                    open(OUT, "| gzip -1c > $path/reads-$haplotype-$fileNumber.fasta.gz");
                }
            }

            close(INP);
        }

        close(OUT);

        open(OUT, "> $path/reads-$haplotype.success");
        close(OUT);
    }

    print STDERR "-- Done!\n";
    print STDERR "--\n";

    #
    #  Figure out what files meryl needs to count
    #

    my @merylInputs;

    foreach my $haplotype (@haplotypes) {
        my @readFiles  = split '\0', $haplotypeReads{$haplotype};
        my $fileNumber = "001";

        while (-e "$path/reads-$haplotype-$fileNumber.fasta.gz") {
            #print STDERR "-- Found input './reads-$haplotype-$fileNumber' in '$path'.\n";

            push @merylInputs, "$haplotype-$fileNumber";
            $fileNumber++;
        }
    }

    #
    #  Emit a script for counting kmers.  A similar version is used in Meryl.pm.
    #

    #  compute kmer size given genome size and error rate
    my $genomeSize = getGlobal("genomeSize");
    my $erate      = 0.001;
    my $merSize    = int(ceil(log($genomeSize * (1 - $erate) / $erate) / log(4)));

    my $mem = getGlobal("merylMemory");
    my $thr = getGlobal("merylThreads");
    my $nJobs;

    open(F, "> $path/meryl-count.sh") or caExit("can't open '$path/meryl-count.sh' for writing: $!", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";

    $nJobs = scalar(@merylInputs);

    for (my $JJ=1; $JJ <= $nJobs; $JJ++) {
        print F "if [ \$jobid -eq $JJ ] ; then\n";
        print F "  batch=\"$merylInputs[$JJ-1]\"\n";
        print F "fi\n";
        print F "\n";
    }

    print F "if [ \$jobid -gt $nJobs ]; then\n";
    print F "  echo Error: Only $nJobs jobs, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e ./reads-\$batch.meryl ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "$bin/meryl k=$merSize threads=$thr memory=$mem \\\n";
    print F "  count \\\n";
    print F "    output ./reads-\$batch.meryl.WORKING \\\n";
    print F "    ./reads-\$batch.fasta.gz \\\n";
    print F "> ./reads-\$batch.meryl.err 2>&1 \\\n";
    print F "&& \\\n";
    print F "mv -f ./reads-\$batch.meryl.WORKING ./reads-\$batch.meryl\n";
    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/meryl-count.sh");
    stashFile("$path/meryl-count.sh");

    #
    #  Emit a script for merging all batches in each haplotype.
    #

    open(F, "> $path/meryl-merge.sh") or caExit("can't open '$path/meryl-count.sh' for writing: $!", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";

    $nJobs = scalar(@haplotypes);

    for (my $JJ=1; $JJ <= $nJobs; $JJ++) {
        print F "if [ \$jobid -eq $JJ ] ; then\n";
        print F "  haplotype=\"$haplotypes[$JJ-1]\"\n";
        print F "fi\n";
        print F "\n";
    }

    print F "$bin/meryl threads=$thr memory=$mem \\\n";
    print F "  union-sum \\\n";
    print F "    output ./reads-\$haplotype.meryl.WORKING \\\n";
    print F "    ./reads-\$haplotype-???.meryl \\\n";
    print F "&& \\\n";
    print F "mv -f ./reads-\$haplotype.meryl.WORKING ./reads-\$haplotype.meryl \\\n";
    print F "&& \\\n";
    print F "rm -rf ./reads-\$haplotype-???.meryl\n";
    print F "\n";
    print F "exit 0\n";

    #
    #  Emit a script for subtracting haplotypes from each other.
    #

    open(F, "> $path/meryl-subtract.sh") or caExit("can't open '$path/meryl-count.sh' for writing: $!", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";

    $nJobs = scalar(@haplotypes);

    for (my $JJ=1; $JJ <= $nJobs; $JJ++) {
        my $oh;

        for (my $KK=1; $KK <= $nJobs; $KK++) {
            if ($JJ != $KK) {
                if (defined($oh)) {
                    $oh .= " ./reads-$haplotypes[$KK-1].meryl";
                } else {
                    $oh  = "./reads-$haplotypes[$KK-1].meryl";
                }
            }
        }

        print F "if [ \$jobid -eq $JJ ] ; then\n";
        print F "  haplotype=\"$haplotypes[$JJ-1]\"\n";
        print F "  otherhaps=\"$oh\"\n";
        print F "fi\n";
        print F "\n";
    }

    print F "$bin/meryl threads=$thr memory=$mem \\\n";
    print F "  difference \\\n";
    print F "    output ./haplotype-\$haplotype.meryl.WORKING \\\n";
    print F "    ./reads-\$haplotype.meryl \\\n";
    print F "    \$otherhaps \\\n";
    print F "&& \\\n";
    print F "mv -f ./haplotype-\$haplotype.meryl.WORKING ./haplotype-\$haplotype.meryl\n";

  finishStage:
    emitStage($asm, "haplotypeCountConfigure");
  allDone:
    stopAfter("meryl");
}



sub haplotypeCountCheck ($) {
    my $asm        = shift @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype/0-kmers";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (skipStage($asm, "haplotype-meryl") == 1);
    goto allDone      if (fileExists("$path/meryl-count.success"));

    #  Scan the script to determine how many jobs there are.

    my @jobs;

    open(F, "< $path/meryl-count.sh") or caExit("can't open '$path/meryl-count.sh' for reading: $!", undef);
    while (<F>) {
        if (m/batch="(.*)"/) {
            push @jobs, $1;
        }
    }
    close(F);

    caExit("failed to find the number of jobs in '$path/meryl-count.sh'", undef)  if (scalar(@jobs) == 0);

    #  Figure out if all the tasks finished correctly.

    my $currentJobID = 1;

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    foreach my $job (@jobs) {
        if      (fileExists("$path/reads-$job.meryl")) {
            push @successJobs, $currentJobID;

        } else {
            $failureMessage .= "--   job $path/meryl-count.sh $currentJobID creating $path/reads-$job.meryl FAILED.\n";
            push @failedJobs, $currentJobID;
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
        emitStage($asm, "haplotype-merylCountCheck", $attempt);

        submitOrRunParallelJob($asm, "meryl", $path, "meryl-count", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " Kmer counting (meryl-count) outputs.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/meryl-count.success") or caExit("can't open '$path/meryl-count.success' for writing: $!", undef);
    close(F);

    stashFile("$path/meryl-count.success");

    generateReport($asm);
    emitStage($asm, "haplotype-merylCountCheck");

  allDone:
    stopAfter("meryl");
}



sub haplotypeMergeCheck ($@) {
    my $asm        = shift @_;
    my @haplotypes =       @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype/0-kmers";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (skipStage($asm, "haplotype-meryl") == 1);
    goto allDone      if (fileExists("$path/meryl-merge.success"));

    #  Figure out if all the tasks finished correctly.  Usually we need to scan the script
    #  to decide how many jobs, but here we just know that there is one job per haplotype.

    my $currentJobID = 1;

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    foreach my $hap (@haplotypes) {
        if      (fileExists("$path/reads-$hap.meryl")) {
            push @successJobs, $currentJobID;

        } else {
            $failureMessage .= "--   job $path/meryl-merge.sh $currentJobID creating $path/reads-$hap.meryl FAILED.\n";
            push @failedJobs, $currentJobID;
        }

        $currentJobID++;
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Kmer counting (meryl-merge) jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Kmer counting (meryl-merge) jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);
        emitStage($asm, "haplotype-merylMergeCheck", $attempt);

        submitOrRunParallelJob($asm, "meryl", $path, "meryl-merge", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " Kmer counting (meryl-merge) outputs.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/meryl-merge.success") or caExit("can't open '$path/meryl-merge.success' for writing: $!", undef);
    close(F);

    stashFile("$path/meryl-merge.success");

    generateReport($asm);
    emitStage($asm, "haplotype-merylMergeCheck");

  allDone:
    stopAfter("meryl");
}



sub haplotypeSubtractCheck ($@) {
    my $asm        = shift @_;
    my @haplotypes =       @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype/0-kmers";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (skipStage($asm, "haplotype-meryl") == 1);
    goto allDone      if (fileExists("$path/meryl-subtract.success"));

    #  Figure out if all the tasks finished correctly.  Usually we need to scan the script
    #  to decide how many jobs, but here we just know that there is one job per haplotype.

    my $currentJobID = 1;

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    foreach my $hap (@haplotypes) {
        if      (fileExists("$path/haplotype-$hap.meryl")) {
            push @successJobs, $currentJobID;

        } else {
            $failureMessage .= "--   job $path/meryl-subtract.sh $currentJobID creating $path/haplotype-$hap.meryl FAILED.\n";
            push @failedJobs, $currentJobID;
        }

        $currentJobID++;
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Kmer counting (meryl-subtract) jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Kmer counting (meryl-subtract) jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);
        emitStage($asm, "haplotype-merylSubtractCheck", $attempt);

        submitOrRunParallelJob($asm, "meryl", $path, "meryl-subtract", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " Kmer counting (meryl-subtract) outputs.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/meryl-subtract.success") or caExit("can't open '$path/meryl-subtract.success' for writing: $!", undef);
    close(F);

    stashFile("$path/meryl-subtract.success");

    generateReport($asm);
    emitStage($asm, "haplotype-merylSubtractCheck");

  allDone:
    stopAfter("meryl");
}



sub haplotypeReadsConfigure ($@) {
    my $asm            = shift @_;
    my $haplotypes     = shift @_;
    my $inputFiles     = shift @_;
    my $bin            = getBinDirectory();
    my $cmd;
    my $path           = "haplotype";

    goto allDone   if (skipStage($asm, "splitHaplotypeConfigure") == 1);
    goto allDone   if (fileExists("$path/splitHaploType.sh"));

    make_path($path)  if (! -d $path);

    #  Canu helpfully prefixes each read with the technology type.  Strip that off.

    my @inputs;

    foreach my $in (@$inputFiles) {
        push @inputs, (split '\0', $in)[1];
    }

    #  splitHaplotype helpfully writes output .fasta.gz where the inpout .meryl files are.
    #  Since this is a quite simple (and the last) step, we'll run it in the root haplotyping/
    #  directory instead of the usual #-stage directory.

    #  Symlink to the meryl databases.

    foreach my $hap (@$haplotypes) {
        if (! -e "haplotype/haplotype-$hap.meryl") {
            symlink("0-kmers/haplotype-$hap.meryl", "haplotype/haplotype-$hap.meryl") or
                caExit("can't make symlink to '0-kmers/haplotype-$hap.meryl' in 'haplotype/haplotype-$hap.meryl': $!", undef);
        }
    }

    #  Now just dump the script.

    my $mem = 48;  #getGlobal("merylMemory");
    my $thr = 16;  #getGlobal("merylThreads");

    open(F, "> $path/splitHaplotype.sh") or caExit("can't open '$path/splitHaplotype.sh' for writing: $!", undef);
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
    print F "if [ -e ./haplotype-unknown.fasta.gz ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "\n";
    print F "$bin/splitHaplotype \\\n";
    print F "  -threads $thr \\\n";
    print F "  -R $_ \\\n"                     foreach (@inputs);       #  One line, yay, but not use of $_.
    print F "  -H ./haplotype-$_.meryl \\\n"   foreach (@$haplotypes);  #  One line, yay, but not use of $_.
    print F "  -A ./haplotype-unknown.fasta.WORKING.gz \\\n";
    print F "> ./splitHaplotype.err 2>&1 \\\n";
    print F "&& \\\n";
    print F "mv -f ./haplotype-unknown.fasta.WORKING.gz ./haplotype-unknown.fasta.gz \\\n";
    print F "&& \\\n";
    print F "exit 0\n";
    close(F);

  finishStage:
    emitStage($asm, "haplotypeReadsConfigure");
  allDone:
    stopAfter("meryl");
}



sub haplotypeReadsCheck ($@) {
    my $asm        = shift @_;
    my @haplotypes =       @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (skipStage($asm, "haplotype-meryl") == 1);
    goto allDone      if (fileExists("$path/haplotyping.success"));

    #  Determine how many jobs we ran.  Usually we need to scan the script to decide how many jobs,
    #  but here we just know that there is exatly one job.  We'll still pretend there are multiple
    #  jobs, just in case someone later implements that.

    my @jobs;

    push @jobs, "1";

    #  Figure out if all the tasks finished correctly.

    my $currentJobID = 1;

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    foreach my $job (@jobs) {
        if      (fileExists("$path/haplotype-unknown.fasta.gz")) {
            push @successJobs, $currentJobID;

        } else {
            $failureMessage .= "--   job $job FAILED.\n";
            push @failedJobs, $currentJobID;
        }

        $currentJobID++;
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Haplotyping jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Haplotyping jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);
        emitStage($asm, "haplotype-merylSubtractCheck", $attempt);

        submitOrRunParallelJob($asm, "meryl", $path, "splitHaplotype", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " haplotyping outputs.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/haplotyping.success") or caExit("can't open '$path/haplotyping.success' for writing: $!", undef);
    close(F);

    stashFile("$path/haplotyping.success");

    generateReport($asm);
    emitStage($asm, "haplotype-merylSubtractCheck");

  allDone:
    stopAfter("meryl");
}
