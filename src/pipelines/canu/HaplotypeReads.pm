
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

package canu::HaplotypeReads;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(haplotypeReadsExist
             haplotypeCountConfigure
             haplotypeCountCheck
             haplotypeMergeCheck
             haplotypeSubtractCheck
             haplotypeReadsConfigure
             haplotypeReadsCheck
             bootstrapHaplotypeAssemblies);

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

    if (scalar(@haplotypes) == 0) {
        return("no-haplotypes");
    }

    #  Check for each file, failing if any one doesn't exist.

    foreach my $haplotype (@haplotypes) {
        return("no")   if (! fileExists("./haplotype/haplotype-$haplotype.fasta.gz"));
    }

    return("no")       if (! fileExists("./haplotype/haplotyping.success"));

    #  Yeah!  Reads exist!

    return("yes");
}



sub haplotypeSplitReads ($$%) {
    my $asm            = shift @_;
    my $merSize        = shift @_;
    my %haplotypeReads =       @_;
    my $bin            = getBinDirectory();
    my $cmd;
    my $path           = "haplotype/0-kmers";

    my @haplotypes     = keys %haplotypeReads;

    #
    #  Check if we're already split.
    #

    my $splitNeeded = 0;

    foreach my $haplotype (@haplotypes) {
        $splitNeeded = 1   if (! fileExists("$path/reads-$haplotype/reads-$haplotype.success"));
    }

    return   if ($splitNeeded == 0);

    #
    #  Blindly split each file into 100 Mbp chunks.
    #    !00x of a 150 Mbp genome ->  150 files.
    #    100x of a   3 Gbp genome -> 3000 files.
    #

    foreach my $haplotype (@haplotypes) {
        my @readFiles     = split '\0', $haplotypeReads{$haplotype};
        my $fileNumber    = "001";
        my $fileLength    = 0;
        my $fileLengthMax = 100000000;

        next  if (fileExists("$path/reads-$haplotype/reads-$haplotype.success"));

        make_path("$path/reads-$haplotype")  if (! -d "$path/reads-$haplotype");

        print STDERR "--\n";
        print STDERR "--  Split haplotype reads for '$haplotype' into multiple files.\n";
        print STDERR "--   --> '$path/reads-$haplotype/reads-$haplotype-$fileNumber.fasta.gz'.\n";
        open(OUT, "| gzip -1c > $path/reads-$haplotype/reads-$haplotype-$fileNumber.fasta.gz");

        foreach my $file (@readFiles) {

            #  If any of the files are links to objects, fetch the object
            #  to local disk and update the name.
            #
            #  A similar blcok is used in SequenceStore.pm and HaplotypeReads.pm (twice).

            if ($file =~ m/dnanexus:(.*)=(.*)/) {
                my $link = $1;
                my $name = $2;

                print STDERR "-- Fetch input file './$name' from object '$link'\n";

                fetchFileFromLink($link, $name);

                $file = "./$name";
            }

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
                    stashFile("$path/reads-$haplotype/reads-$haplotype-$fileNumber.fasta.gz");

                    $fileNumber++;
                    $fileLength = 0;

                    print STDERR "--   --> '$path/reads-$haplotype/reads-$haplotype-$fileNumber.fasta.gz'.\n";
                    open(OUT, "| gzip -1c > $path/reads-$haplotype/reads-$haplotype-$fileNumber.fasta.gz");
                }
            }

            close(INP);
        }

        close(OUT);

        stashFile("$path/reads-$haplotype/reads-$haplotype-$fileNumber.fasta.gz");

        open(OUT, "> $path/reads-$haplotype/reads-$haplotype.success");
        print OUT "$fileLengthMax\n";    #  Not really used
        print OUT "$fileNumber\n";       #
        close(OUT);

        stashFile("$path/reads-$haplotype/reads-$haplotype.success");
    }
}



sub haplotypeCountConfigure ($%) {
    my $asm            = shift @_;
    my %haplotypeReads =       @_;
    my $bin            = getBinDirectory();
    my $cmd;
    my $path           = "haplotype/0-kmers";

    my @haplotypes     = sort keys %haplotypeReads;

    goto allDone   if (fileExists("$path/meryl-count.sh") &&
                       fileExists("$path/meryl-merge.sh") &&
                       fileExists("$path/meryl-subtract.sh"));

    make_path($path)  if (! -d $path);

    #
    #  Pick an appropriate mer size based on genome size.
    #

    my $genomeSize = getGlobal("genomeSize");
    my $erate      = 0.001;
    my $merSize    = int(ceil(log($genomeSize * (1 - $erate) / $erate) / log(4)));


    #
    #  Split reads if we need to.
    #

    haplotypeSplitReads($asm, $merSize, %haplotypeReads);

    #
    #  Figure out what files meryl needs to count
    #

    my %merylInputs;
    my %merylFiles;
    my $totalFiles;

    print STDERR "--\n";

    foreach my $haplotype (@haplotypes) {
        my @readFiles  = split '\0', $haplotypeReads{$haplotype};
        my $fileNumber = "001";

        while (fileExists("$path/reads-$haplotype/reads-$haplotype-$fileNumber.fasta.gz")) {
            $merylInputs{$haplotype} .= "$haplotype-$fileNumber\0";
            $merylFiles{$haplotype}++;

            $totalFiles++;

            $fileNumber++;
        }

        print STDERR "-- haplotype $haplotype has $merylFiles{$haplotype} input files.\n";
    }

    #
    #  Figure out batch sizes for meryl.   - the number of jobs per haplotype.
    #
    #  -- Up to  5 files/job with  2 jobs ->   10 files.
    #  -- Up to  8 files/job with  4 jobs ->   32 files.
    #  -- Up to 11 files/job with  6 jobs ->   66 files.
    #  -- Up to 14 files/job with  8 jobs ->  112 files.
    #  -- Up to 17 files/job with 10 jobs ->  170 files.
    #  -- Up to 20 files/job with 12 jobs ->  240 files.
    #  -- Up to 23 files/job with 14 jobs ->  322 files.
    #  -- Up to 26 files/job with 16 jobs ->  416 files.
    #  -- Up to 29 files/job with 18 jobs ->  522 files.
    #  -- Up to 32 files/job with 20 jobs ->  640 files.
    #  -- Up to 35 files/job with 22 jobs ->  770 files.
    #  -- Up to 38 files/job with 24 jobs ->  912 files.
    #  -- Up to 41 files/job with 26 jobs -> 1066 files.
    #  -- Up to 44 files/job with 28 jobs -> 1232 files.
    #  -- Up to 47 files/job with 30 jobs -> 1410 files.
    #  -- Up to 50 files/job with 32 jobs -> 1600 files.
    #  -- Up to 53 files/job with 34 jobs -> 1802 files.
    #  -- Up to 56 files/job with 36 jobs -> 2016 files.
    #  -- Up to 59 files/job with 38 jobs -> 2242 files.
    #  -- Up to 62 files/job with 40 jobs -> 2480 files.
    #  -- Up to 65 files/job with 42 jobs -> 2730 files.
    #  -- Up to 68 files/job with 44 jobs -> 2992 files.
    #  -- Up to 71 files/job with 46 jobs -> 3266 files.
    #  -- Up to 74 files/job with 48 jobs -> 3552 files.
    #  -- Up to 77 files/job with 50 jobs -> 3850 files.
    #  -- Up to 80 files/job with 52 jobs -> 4160 files.
    #  -- Up to 83 files/job with 54 jobs -> 4482 files.
    #  -- Up to 86 files/job with 56 jobs -> 4816 files.
    #  -- Up to 89 files/job with 58 jobs -> 5162 files.
    #  -- Up to 92 files/job with 60 jobs -> 5520 files.
    #  -- Up to 95 files/job with 62 jobs -> 5890 files.
    #  -- Up to 98 files/job with 64 jobs -> 6272 files.
    #

    my $nFilesPerJob = 5;
    my $nJobsMax     = scalar(@haplotypes);

    #print STDERR "-- Up to $nFilesPerJob files/job with $nJobsMax jobs -> ", $nFilesPerJob * $nJobsMax, " files.\n";

    while ($nFilesPerJob * $nJobsMax < $totalFiles) {
        $nFilesPerJob += 3;
        $nJobsMax     += scalar(@haplotypes);

        #print STDERR "-- Up to $nFilesPerJob files/job with $nJobsMax jobs -> ", $nFilesPerJob * $nJobsMax, " files.\n";
    }

    print STDERR "-- Will use $nJobsMax jobs with up to $nFilesPerJob files each.\n";

    #
    #  Assign read files to jobs.  Divide the files for each haplotype into groups
    #  of some size, taking care to not mix up files between haplotypes!
    #

    my @merylInputs;             #  A list of the inputs to each counting job.
    my @merylOutputs;            #  The output of each counting job.
    my %merylJobsPerHaplotype;   #  The number of counting jobs for each haplotype.

    my $jid = "001";

    foreach my $haplotype (@haplotypes) {
        my @files  = split '\0', $merylInputs{$haplotype};
        my $nFiles =             $merylFiles{$haplotype};

        while (scalar(@files) > 0) {
            my $files = undef;

            for (my $n=0; $n<$nFilesPerJob; $n++) {
                last   if (!defined($files[0]));

                if (! defined($files)) {
                    $files  = "  batch=\"./reads-$haplotype/reads-$files[0].fasta.gz \\\n";
                } else {
                    $files .= "         ./reads-$haplotype/reads-$files[0].fasta.gz \\\n";
                }

                shift @files;
            }

            $files =~ s/gz\s\\\n$/gz"/;   #  Remove the " \" at the end of the string, and add a terminating quote.

            push @merylInputs, $files;               #  Save the list of inputs.
            push @merylOutputs, "$haplotype-$jid";   #  Save the output name.

            if (defined($merylJobsPerHaplotype{$haplotype})) {
                $merylJobsPerHaplotype{$haplotype} .= " $jid";
            } else {
                $merylJobsPerHaplotype{$haplotype}  = "$jid";
            }

            $jid++;
        }
    }

    #
    #  Memory requirements are dependent only on $nFilesPerJob...unless you
    #  change the size of each file, so don't do that!  Well, OK, and merSize.
    #
    #  Observed memory usage for 19-mers on a 150 Gbp genome:
    #     1 -> 1.0 GB
    #     9 -> 3.4 GB
    #    15 ->
    #    23 -> 8.2 GB
    #    28 -> 8.0 GB
    #
    #  So roughly 300 MB per 100 MB input.  We'll overestimate the memory we
    #  request, and tell meryl the smaller size.  So if we do screw up the
    #  300 MB per 100 Mb claim, meryl itself will dump partial results and
    #  merge at the end.
    #

    my $thr = getGlobal("merylThreads");
    my $mem = 2 + ceil(0.333 * $nFilesPerJob + 0.5);

    open(F, "> $path/meryl-count.memory") or caExit("can't open '$path/meryl-count.memory' for writing: $!", undef);
    print F "$mem\n";
    close(F);
    stashFile("$path/meryl-count.memory");

    $mem = 1.0 + 0.333 * $nFilesPerJob;

    #
    #  Make a script for counting.  All the counting jobs are the same, no need to
    #  distinguish between haplotypes here.
    #

    setGlobal("merylMemory", $mem);

    open(F, "> $path/meryl-count.sh") or caExit("can't open '$path/meryl-count.sh' for writing: $!", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";

    my $nJobs = scalar(@merylInputs);

    for (my $JJ=1; $JJ <= $nJobs; $JJ++) {
        print F "if [ \$jobid -eq $JJ ] ; then\n";
        print F "$merylInputs[$JJ-1]\n";
        print F "  output=\"$merylOutputs[$JJ-1]\"\n";
        print F "fi\n";
        print F "\n";
    }

    print F "\n";
    print F "if [ \$jobid -gt $nJobs ]; then\n";
    print F "  echo Error: Only $nJobs jobs, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";

    print F "\n";
    print F "#  If the meryl database exists, we're done.\n";
    print F "\n";
    print F "if [ -e ./reads-\$output.meryl/merylIndex ] ; then\n";
    print F "  echo Kmers for batch \$output exist.\n";
    print F "  exit 0\n";
    print F "fi\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  If the meryl output exists in the object store, we're also done.\n";
        print F "\n";
        print F fileExistsShellCode("exist1", "$path", "reads-\$output.meryl.tar.gz");
        print F "if [ \$exist1 = true ] ; then\n";
        print F "  echo Kmers for batch \$output exist in the object store.\n";
        print F "  exit 0\n";
        print F "fi\n";
    }

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Nope, not done.  Fetch input files.\n";
        print F "\n";
        print F "for file in \$batch ; do\n";
        print F fetchFileShellCode($path, "\$file", "  ");
        print F "done\n";
    }

    print F "\n";
    print F "#  And compute.\n";
    print F "\n";
    print F "\$bin/meryl k=$merSize threads=$thr memory=$mem \\\n";
    print F "  count \\\n";
    print F "    output ./reads-\$output.meryl.WORKING \\\n";
    print F "    \$batch \\\n";
    print F "&& \\\n";
    print F "mv -f ./reads-\$output.meryl.WORKING ./reads-\$output.meryl\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Save those precious results!\n";
        print F "\n";
        print F stashMerylShellCode($path, "reads-\$output.meryl", "  ");
    }

    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/meryl-count.sh");
    stashFile("$path/meryl-count.sh");

    #
    #  Emit a script for merging all batches in each haplotype.  Unlike counting,
    #  we definitely need to take the haplotype the job is for into consideration.
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

    for (my $JJ=1; $JJ <= $nJobs; $JJ++) {           #  Set the name of the haplotype this job
        my $hap = $haplotypes[$JJ-1];                #  is merging kmers for.

        print F "if [ \$jobid -eq $JJ ] ; then\n";
        print F "  haplotype=\"$hap\"\n";
        print F "  batches=\"$merylJobsPerHaplotype{$hap}\"\n";
        print F "fi\n";
        print F "\n";
    }

    print F "\n";
    print F "#  If the meryl database exists, we're done.\n";
    print F "\n";
    print F "if [ -e ./reads-\$haplotype.meryl/merylindex ] ; then\n";
    print F "  echo Merged kmers for haplotype \$haplotype exist.\n";
    print F "  exit 0\n";
    print F "fi\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  If the meryl output exists in the object store, we're also done.\n";
        print F "\n";
        print F fileExistsShellCode("exist1", "$path", "reads-\$haplotype.meryl.tar.gz");
        print F "if [ \$exist1 = true ] ; then\n";
        print F "  echo Merged kmers for haplotype \$haplotype exist in the object store.\n";
        print F "  exit 0\n";
        print F "fi\n";
    }

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Nope, not done.  Fetch the input databases.\n";
        print F "\n";
        print F "for batch in \$batches ; do\n";
        print F fetchMerylShellCode($path, "reads-\$haplotype-\$batch.meryl", "  ");
        print F "done\n";
    }

    print F "\n";
    print F "#  And compute.\n";
    print F "\n";
    print F "\$bin/meryl threads=$thr memory=$mem \\\n";
    print F "  union-sum \\\n";
    print F "    output ./reads-\$haplotype.meryl.WORKING \\\n";
    print F "    ./reads-\$haplotype-???.meryl \\\n";
    print F "&& \\\n";
    print F "mv -f ./reads-\$haplotype.meryl.WORKING ./reads-\$haplotype.meryl \\\n";
    print F "&& \\\n";
    print F "rm -rf ./reads-\$haplotype-???.meryl\n";
    print F "\n";
    print F "\$bin/meryl threads=$thr memory=$mem \\\n";
    print F "  statistics \\\n";
    print F "    ./reads-\$haplotype.meryl \\\n";
    print F "> ./reads-\$haplotype.statistics \\\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Save those precious results!\n";
        print F "\n";
        print F stashFileShellCode($path, "reads-\$haplotype.statistics", "");
        print F stashMerylShellCode($path, "reads-\$haplotype.meryl", "");
    }

    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/meryl-merge.sh");
    stashFile("$path/meryl-merge.sh");

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

    print F "\n";
    print F "#  If the meryl database exists, we're done.\n";
    print F "\n";
    print F "if [ -e ./haplotype-\$haplotype.meryl ] ; then\n";
    print F "  echo Kmers for haplotype \$haplotype exist.\n";
    print F "  exit 0\n";
    print F "fi\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  If the meryl output exists in the object store, we're also done.\n";
        print F "\n";
        print F fileExistsShellCode("exist1", "$path", "haplotype-\$haplotype.meryl.tar.gz");
        print F "if [ \$exist1 = true ] ; then\n";
        print F "  echo Kmers for haplotype \$haplotype exist in the object store.\n";
        print F "  exit 0\n";
        print F "fi\n";
    }

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Fetch _all_ the haplotype kmers.\n";
        print F "\n";
        print F fetchMerylShellCode($path, "reads-$_.meryl", "")   foreach (@haplotypes);   #  One line, yay, but not use of $_.
    }

    print F "\n";
    print F "#  Subtract all the other haplotypes from ours.\n";
    print F "\n";
    print F "\$bin/meryl threads=$thr memory=$mem \\\n";
    print F "  difference \\\n";
    print F "    output ./haplotype-\$haplotype.meryl.WORKING \\\n";
    print F "    ./reads-\$haplotype.meryl \\\n";
    print F "    \$otherhaps \\\n";
    print F "&& \\\n";
    print F "mv -f ./haplotype-\$haplotype.meryl.WORKING ./haplotype-\$haplotype.meryl\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Save those precious results!\n";
        print F "\n";
        print F stashMerylShellCode($path, "haplotype-\$haplotype.meryl", "");
    }

    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/meryl-subtract.sh");
    stashFile("$path/meryl-subtract.sh");

  finishStage:
    resetIteration("haplotypeCountConfigure");

  allDone:
    stopAfter("meryl-configure");
}



sub haplotypeCountCheck ($) {
    my $asm        = shift @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype/0-kmers";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (fileExists("$path/meryl-count.success"));

    #  Scan the script to determine how many jobs there are.

    fetchFile("$path/meryl-count.sh");

    my @jobs;

    open(F, "< $path/meryl-count.sh") or caExit("can't open '$path/meryl-count.sh' for reading: $!", undef);
    while (<F>) {
        if (m/output="(.*-\d+)"$/) {
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
        if ((fileExists("$path/reads-$job.meryl")) ||
            (fileExists("$path/reads-$job.meryl.tar.gz"))) {
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

        #  One off hack for setting memory here

        fetchFile("$path/meryl-count.memory");

        my $mem = 0;
        open(F, "< $path/meryl-count.memory") or caExit("can't open '$path/meryl-count.memory' for reading: $!", undef);
        $mem = <F>;  chomp $mem;
        close(F);

        setGlobal("merylMemory", $mem);

        #  And _now_ run some jobs.

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
    resetIteration("haplotype-merylCountCheck");

  allDone:
    stopAfter("meryl-count");
}



sub haplotypeMergeCheck ($@) {
    my $asm        = shift @_;
    my @haplotypes =       @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype/0-kmers";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (fileExists("$path/meryl-merge.success"));

    fetchFile("$path/meryl-merge.sh");

    #  Figure out if all the tasks finished correctly.  Usually we need to scan the script
    #  to decide how many jobs, but here we just know that there is one job per haplotype.

    my $currentJobID = 1;

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    foreach my $hap (@haplotypes) {
        if ((fileExists("$path/reads-$hap.meryl")) ||
            (fileExists("$path/reads-$hap.meryl.tar.gz"))) {
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
    resetIteration("haplotype-merylMergeCheck");

  allDone:
    stopAfter("meryl-merge");
}



sub haplotypeSubtractCheck ($@) {
    my $asm        = shift @_;
    my @haplotypes =       @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype/0-kmers";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (fileExists("$path/meryl-subtract.success"));

    fetchFile("$path/meryl-subtract.sh");

    #  Figure out if all the tasks finished correctly.  Usually we need to scan the script
    #  to decide how many jobs, but here we just know that there is one job per haplotype.

    my $currentJobID = 1;

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    foreach my $hap (@haplotypes) {
        if ((fileExists("$path/haplotype-$hap.meryl")) ||
            (fileExists("$path/haplotype-$hap.meryl.tar.gz"))) {
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
    resetIteration("haplotype-merylSubtractCheck");

    #  Remove databases.

    print STDERR "--\n";

    foreach my $hap (@haplotypes) {
        if (getGlobal("saveMerCounts") == 0) {
            print STDERR "-- Removing meryl database '$path/reads-$hap.meryl'.\n";
            remove_tree("$path/reads-$hap.meryl")
        } else {
            print STDERR "-- Meryl database '$path/reads-$hap.meryl' saved because 'saveMerCounts=true'.\n";
        }
    }

  allDone:
    stopAfter("meryl-subtract");
}



sub haplotypeReadsConfigure ($@) {
    my $asm            = shift @_;
    my $haplotypes     = shift @_;
    my $inputFiles     = shift @_;
    my $bin            = getBinDirectory();
    my $cmd;
    my $path           = "haplotype";

    goto allDone   if (fileExists("$path/splitHaplotype.sh"));

    make_path($path)  if (! -d $path);

    #  Canu helpfully prefixes each read with the technology type.  Strip that off.

    my @inputs;

    foreach my $in (@$inputFiles) {
        push @inputs, (split '\0', $in)[1];
    }

    #  Dump the script.

    my $mem = getGlobal("hapMemory");
    my $thr = getGlobal("hapThreads");

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
    print F "#  If the unknown haplotype assignment exists, we're done.\n";
    print F "\n";
    print F "if [ -e ./haplotype.log ] ; then\n";
    print F "  echo Read to haplotype assignment already exists.\n";
    print F "  exit 0\n";
    print F "fi\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  If the unknown haplotype assignment exists in the object store, we're also done.\n";
        print F "\n";
        print F fileExistsShellCode("exist1", "$path", "haplotype.log");
        print F "if [ \$exist1 = true ] ; then\n";
        print F "  echo Read to haplotype assignment exists in the object store.\n";
        print F "  exit 0\n";
        print F "fi\n";
    }

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Fetch _all_ the haplotype kmers.  These need to be in 0-kmers.\n";
        print F "\n";
        print F "mkdir 0-kmers\n";
        print F "cd    0-kmers\n";
        print F "\n";
        print F fetchMerylShellCode("$path/0-kmers", "haplotype-$_.meryl", "")   foreach (@$haplotypes);  #  One line, yay, but not use of $_.
        print F "\n";
        print F "cd ..\n";
    }

    #  Fetch the reads from the object store.
    #  A similar blcok is used in SequenceStore.pm and HaplotypeReads.pm (twice).

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#\n";
        print F "#  Fetch all the unlcassified input reads.\n";
        print F "#\n";

        for (my $ff=0; $ff < scalar(@inputs); $ff++) {
            if ($inputs[$ff] =~ m/dnanexus:(.*)=(.*)/) {
                my $link = $1;
                my $name = $2;

                print F fetchFileFromLinkShellCode($link, $name, "");

                $inputs[$ff] = $name;
            }
        }

        print F "\n";
    }

    my $minReadLength = getGlobal("minReadLength");

    print F "\n";
    print F "#  Assign reads to haplotypes.\n";
    print F "\n";
    print F "\$bin/splitHaplotype \\\n";
    print F "  -cl $minReadLength \\\n";
    print F "  -threads $thr \\\n";
    print F "  -memory  $mem \\\n";
    print F "  -R $_ \\\n"                                                                                   foreach (@inputs);
    print F "  -H ./0-kmers/haplotype-$_.meryl ./0-kmers/reads-$_.statistics ./haplotype-$_.fasta.gz \\\n"   foreach (@$haplotypes);
    print F "  -A ./haplotype-unknown.fasta.gz \\\n";
    print F "> haplotype.log.WORKING \\\n";
    print F "&& \\\n";
    print F "mv -f ./haplotype.log.WORKING ./haplotype.log\n";

    if (defined(getGlobal("objectStore"))) {
        print F "\n";
        print F "#  Save those precious results!\n";
        print F "\n";
        print F stashFileShellCode($path, "haplotype-$_.fasta.gz", "")   foreach (@$haplotypes);
        print F stashFileShellCode($path, "haplotype-unknown.fasta.gz", "");
        print F stashFileShellCode($path, "haplotype.log", "");
    }

    print F "\n";
    print F "exit 0\n";
    close(F);

    makeExecutable("$path/splitHaplotype.sh");
    stashFile("$path/splitHaplotype.sh");

  finishStage:
    resetIteration("haplotypeReadsConfigure");

  allDone:
    stopAfter("haplotype-configure");
}



#  Checks for presence of haplotype.log, and if that exists,
#  declares the haplotyping stage done and creates haplotping.success.
#
#  If haplotyping failed, the (partial) log will be in haplotype.log.WORKING.
#
sub haplotypeReadsCheck ($@) {
    my $asm        = shift @_;
    my @haplotypes =       @_;
    my $attempt    = getGlobal("canuIteration");

    my $path       = "haplotype";

    my $bin        = getBinDirectory();
    my $cmd;

    #  Check if we're known to be done.

    goto allDone      if (fileExists("$path/haplotyping.success"));

    fetchFile("$path/splitHaplotype.sh");

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
        if      (fileExists("$path/haplotype.log")) {
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
            caExit(undef, "$path/haplotype.log.WORKING");
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Haplotyping jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "hap", $path, "splitHaplotype", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " haplotyping outputs.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/haplotyping.success") or caExit("can't open '$path/haplotyping.success' for writing: $!", undef);
    close(F);

    stashFile("$path/haplotyping.success");

    generateReport($asm);
    resetIteration("haplotype-reads");

    #remove_tree("$path/$name")   if (getGlobal("saveMerCounts") == 0);

  allDone:
    stopAfter("haplotype");
}



sub bootstrapHaplotypeAssemblies ($@) {
    my $bin        = getBinDirectory();
    my $asm        = shift @_;
    my @haplotypes =       @_;

    my $techtype   = removeHaplotypeOptions();
    my @options    = getCommandLineOptions();

    #  Find the maximum length of haplotype names, to make the output pretty.

    my $displLen = 0;

    foreach my $haplotype (@haplotypes) {
        my $hapLen = length($haplotype);
        $displLen = ($displLen < $hapLen) ? $hapLen : $displLen;
    }

    #  Decide if we should use or ignore the unassigned reads, and if we should
    #  even bother assembling.

    fetchFile("");

    my %hapReads;
    my %hapBases;

    my $totReads = 0;
    my $totBases = 0;

    open(F, "< haplotype/haplotype.log") or caExit("can't open 'haplotype/haplotype.log' for reading: $!", undef);
    while (<F>) {
        if (m/(\d+)\s+reads\s+(\d+)\s+bases\s+written\s+to\s+haplotype\s+file\s+.*haplotype-(\w+).fasta.gz/) {
            $hapReads{$3} = $1;
            $hapBases{$3} = $2;

            $totReads += $1   if ($3 ne "unknown");
            $totBases += $2   if ($3 ne "unknown");
        }
        if (m/(\d+)\s+reads\s+(\d+)\s+bases\s+filtered\s+for\s+being\s+too\s+short/) {
            $hapReads{"short"} = $1;
            $hapBases{"short"} = $2;
        }
    }
    close(F);

    print STDERR "--\n";
    foreach my $haplotype (@haplotypes) {
        printf STDERR "-- Found   %8d reads and %12d bases for haplotype $haplotype.\n", $hapReads{$haplotype}, $hapBases{$haplotype};
    }
    printf STDERR "-- Found   %8d reads and %12d bases assigned to no haplotype.\n", $hapReads{"unknown"}, $hapBases{"unknown"};
    printf STDERR "-- Ignored %8d reads and %12d bases because they were short.\n",  $hapReads{"short"},   $hapBases{"short"};

    #  Ignore the unknown reads if there aren't that many.

    my $unknownFraction =  getGlobal("hapUnknownFraction");
    my $withUnknown = (($totBases > 0) && ($hapBases{"unknown"} / $totBases < $unknownFraction)) ? 0 : 1;

    if ($withUnknown == 0) {
        print STDERR "--\n";
        print STDERR "-- Fewer than " . $unknownFraction*100 . " % of bases in unassigned reads; don't use them in assemblies.\n";
    } else {
        print STDERR "--\n";
        print STDERR "-- More than " .  $unknownFraction*100 . " % of bases in unassigned reads; including them in assemblies.\n";
    }

    #  For each haplotype, emit a script to run the assembly.

    print STDERR "--\n";
    print STDERR "-- Haplotype assembly commands:\n";

    foreach my $haplotype (@haplotypes) {
        my $hs = substr("$haplotype" . " " x $displLen, 0, $displLen);

        print STDERR "--   ./$asm-haplotype$haplotype.sh\n";
        #print STDERR "--   $rootdir/$asm-haplotype$haplotype.sh\n";

        open(F, "> ./$asm-haplotype$haplotype.sh");
        print F "#!/bin/sh\n";
        print F "\n";

        if (defined(getGlobal("objectStore"))) {
            print F "\n";
            print F "#  Fetch the haplotyped reads.  This is just a bit weird.\n";
            print F "#  The fetch (boilerplate) only works from within a subdirectory,\n";
            print F "#  so we must cd into it first, fetch, the go back to the root.\n";
            print F "\n";
            print F "mkdir -p haplotype\n";
            print F "cd       haplotype\n";
            print F fetchFileShellCode("haplotype", "haplotype-$haplotype.fasta.gz", "");
            print F fetchFileShellCode("haplotype", "haplotype-unnown.fasta.gz", "")      if ($withUnknown);
            print F "cd ..\n";
        }

        print F "\n";
        print F "$bin/canu \\\n";
        print F "  -p $asm-haplotype$haplotype \\\n";
        print F "  -d $asm-haplotype$haplotype \\\n";
        print F "  $_ \\\n"   foreach (@options);
        print F "  $techtype ./haplotype/haplotype-$haplotype.fasta.gz \\\n";
        print F "  $techtype ./haplotype/haplotype-unknown.fasta.gz \\\n"     if ($withUnknown);
        print F "> ./$asm-haplotype$haplotype.out 2>&1\n";
        print F "\n";
        print F "exit 0\n";
        print F "\n";
        close(F);

        makeExecutable("./$asm-haplotype$haplotype.sh");
    }

    #  Fail if too many unassigned reads.

    if ($totBases == 0) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  No reads assigned to haplotypes.  Assemblies not started.\n";
        print STDERR "-- ERROR:\n";
    }

    elsif ($hapBases{"unknown"} / $totBases > 0.50) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Too many bases in unassigned reads.  Assemblies not started.\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  If you run them manually, note that the unassigned reads\n";
        print STDERR "-- ERROR:  are included in ALL assemblies.\n";
        print STDERR "-- ERROR:\n";
    }

    #  Finished.  Let the caller submit the assemblies or not.
}
