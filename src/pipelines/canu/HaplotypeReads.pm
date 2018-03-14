
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
 #    Brian P. Walenz beginning on 2018-FEB-08
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::HaplotypeReads;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(checkHaplotypeReads getHaplotypes haplotypeConfigure haplotypeCheck dumpHaplotypeReads);

use strict;

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Configure;
use canu::Execution;
use canu::Gatekeeper;
use canu::Report;
use canu::Grid_Cloud;

sub getHaplotypes($) {
   my @haplotypes;
   my $base = shift @_;

   opendir(DIR, $base);
   while (my $file = readdir(DIR)) {
      next unless (-d "$base/$file");
      next unless ($file =~ m/0-mer/);
      if ($file =~ m/0-mercounts-(\S*)/) {
         push @haplotypes, $1;
      }
   }
   closedir(DOR);

   return @haplotypes;
}

sub checkHaplotypeReads($$) {
    my $asm  = shift @_;
    my $base = shift @_;

    my @haplotypes = getHaplotypes("haplotype");
    my $allSequencesDone = 1;

    push @haplotypes, "unknown";
    foreach my $haplotype (@haplotypes) {
       if (!sequenceFileExists("$asm.${haplotype}Reads.fasta.gz")) {
          $allSequencesDone = 0;
       }
    }

    return $allSequencesDone;
}

#  Return the number of correction jobs.
#
sub computeNumberOfHaplotypeJobs ($) {
    my $asm     = shift @_;
    my $nJobs   = 0;
    my $nPerJob = 0;

    my $nPart    = getGlobal("corPartitions");
    my $nReads   = getNumberOfReadsInStore("hap", $asm);

    caExit("didn't find any reads in store 'haplotype/$asm.gkpStore'?", undef)  if ($nReads == 0);

    $nPerJob     = int($nReads / $nPart + 1);
    $nPerJob     = getGlobal("corPartitionMin")  if ($nPerJob < getGlobal("corPartitionMin"));

    for (my $j=1; $j<=$nReads; $j += $nPerJob) {  #  We could just divide, except for rounding issues....
        $nJobs++;
    }

    return($nJobs, $nPerJob);
}

sub estimateMemoryNeededForHaplotypeJobs ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $tag     = "hap";

    my $path    = "haplotype/2-correction";

    my $memEst   = 0;

    return   if (defined(getGlobal("corMemory")));

    my @haplotypes = getHaplotypes("haplotype");
    my $merSize    = getGlobal("${tag}OvlMerSize");

    foreach my $haplotype (@haplotypes) {
       fetchFile("haplotype/0-mercounts-$haplotype/$haplotype.ms$merSize.only.mcdat");

       if (-e "haplotype/0-mercounts-$haplotype/$haplotype.ms$merSize.only.mcdat") {
          my $size = -s "haplotype/0-mercounts-$haplotype/$haplotype.ms$merSize.only.mcdat";
          $memEst = int($size / 1073741824.0 + 0.5) * 2;
        }
        close(F);
    }

    if ($memEst == 0) {
        $memEst = 12;
    }

    setGlobal("corMemory", $memEst);

    my $err;
    my $all;

    ($err, $all) = getAllowedResources("", "cor", $err, $all, 0);

    print STDERR "--\n";
    print STDERR $all;
    print STDERR "--\n";
}

sub haplotypeConfigure ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "haplotype";
    my $path    = "haplotype/2-haplotype";
    my $tag     = "hap";

    make_path("$path")  if (! -d "$path");

    goto allDone   if (skipStage($asm, "hap-haplotypeConfigure") == 1);
    goto allDone   if (fileExists("$path/haplotypeReads.sh"));              #  Jobs created

    fetchStore("./correction/$asm.gkpStore");

    estimateMemoryNeededForHaplotypeJobs($asm);

    my ($nJobs, $nPerJob)  = computeNumberOfHaplotypeJobs($asm);  #  Does math based on number of reads and parameters.

    my $nReads             = getNumberOfReadsInStore("hap", $asm);

    open(F, "> $path/haplotypeReads.sh") or caExit("can't open '$path/haplotypeReads.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("$base/$asm.gkpStore", "$base/2-correction", "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    print F "if [ \$jobid -gt $nJobs ]; then\n";
    print F "  echo Error: Only $nJobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";

    my  $bgnID   = 1;
    my  $endID   = $bgnID + $nPerJob - 1;
    my  $jobID   = 1;

    while ($bgnID < $nReads) {
        $endID  = $bgnID + $nPerJob - 1;
        $endID  = $nReads  if ($endID > $nReads);

        print F "if [ \$jobid -eq $jobID ] ; then\n";
        print F "  bgn=$bgnID\n";
        print F "  end=$endID\n";
        print F "fi\n";

        $bgnID = $endID + 1;
        $jobID++;
    }

    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ -e \"./results/\$jobid.success\" ] ; then\n";
    print F "  echo Job finished successfully.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d \"./results\" ] ; then\n";
    print F "  mkdir -p \"./results\"\n";
    print F "fi\n";
    print F "\n";

    print F fetchStoreShellCode("haplotype/$asm.gkpStore", "haplotype/2-haplotype", "");
    print F "\n";
    # fetch the filter files

    print F "gkpStore=\"../$asm.gkpStore\"\n";
    print F "\n";

    my $stageDir = getGlobal("stageDirectory");

    if (defined($stageDir)) {
        print F "if [ ! -d $stageDir ] ; then\n";
        print F "  mkdir -p $stageDir\n";
        print F "fi\n";
        print F "\n";
        print F "mkdir -p $stageDir/$asm.gkpStore\n";
        print F "\n";
        print F "echo Start copy at `date`\n";
        print F "cp -p \$gkpStore/info      $stageDir/$asm.gkpStore/info\n";
        print F "cp -p \$gkpStore/libraries $stageDir/$asm.gkpStore/libraries\n";
        print F "cp -p \$gkpStore/reads     $stageDir/$asm.gkpStore/reads\n";
        print F "cp -p \$gkpStore/blobs.*   $stageDir/$asm.gkpStore/\n";
        print F "echo Finished   at `date`\n";
        print F "\n";
        print F "gkpStore=\"$stageDir/$asm.gkpStore\"\n";
        print F "\n";
    }

    print F "\n";
    print F "\$bin/gatekeeperDumpFASTQ \\\n";
    print F "  -G \$gkpStore \\\n";
    print F "  -o ./results/\$jobid.fasta.WORKING \\\n";
    print F "  -r \$bgn-\$end \\\n";
    print F "  -fasta -nolibname -noreadname \\\n";
    print F "  > ./results/\$jobid.err 2>&1 \\\n";
    print F "&& \\\n";
    print F "mv ./results/\$jobid.fasta.WORKING.fasta ./results/\$jobid.fasta \\\n";
    print F "\n";

    my @haplotypes = getHaplotypes($base);
    my $merSize    = getGlobal("${tag}OvlMerSize");
    my $lo         = 0;
    my $hi         = 1000;

    foreach my $haplotype (@haplotypes) {
       fetchFile("$base/0-mercounts-$haplotype/$haplotype.ms$merSize.threshold");
       open(T, "< haplotype/0-mercounts-$haplotype/$haplotype.ms$merSize.threshold") or caExit("can't open haplotype/0-mercounts-$haplotype/$haplotype.ms$merSize.threshold", undef);
       my $rn = <T>;
       ($rn =~ m/^(\d+)\s+(\d+)$/);
       $lo = $1;
       $hi = $2;
       close(T);

       print F "\$bin/simple-dump \\\n";
       print F "  -m $merSize \\\n";
       print F "  -s ../0-mercounts-$haplotype/$haplotype.ms$merSize.only \\\n";
       print F "  -l $lo -h $hi -f ./results/\$jobid.fasta \\\n";
       print F "  > ./results/\$jobid.$haplotype.WORKING \\\n";
       print F "&& \\\n";
       print F "mv ./results/\$jobid.$haplotype.WORKING ./results/\$jobid.$haplotype \\\n";
       print F "\n";
    }
    print F "join results/\$jobid.haplotype* |awk '{print \$1\" \"\$4\" \"\$5\" \"\$8\" \"\$NF}'|awk '{print \$1\" \"\$3/\$2\" \"\$NF/\$(NF-1)}'|awk '{if (\$2 == 0 && \$3 == 0) print \$1\" ambiguous\"; else if (\$2 < \$3) { print \$1\" haplotype2\"; } else if (\$2 > \$3) { print \$1\" haplotype1\"; } else print \$1\" ambiguous\"}'|awk '{print \$NF}'|sort |uniq -c\n";

    # now classify and dump the fasta
    print F "\$bin/splitHaplotype \\\n";
    print F "  -G \$gkpStore \\\n";
    print F "  -p results/\$jobid \\\n";
    print F "  -h " . join(" ", @haplotypes) . "\\\n";
    print F "  -cr 1 -cl " . getGlobal("minReadLength") . " \\\n";
    print F "  -b \$bgn -e \$end \\\n";
    print F "&& \\\n";
    print F "touch ./results/\$jobid.success \\\n";
    print F "\n";

    if (defined($stageDir)) {
        print F "rm -rf $stageDir/$asm.gkpStore\n";   #  Prevent accidents of 'rm -rf /' if stageDir = "/".
        print F "rmdir  $stageDir\n";
        print F "\n";
    }

    print F stashFileShellCode("$path", "results/\$jobid.cns", "");

    print F "\n";
    print F "exit 0\n";

    close(F);

    makeExecutable("$path/haplotypeReads.sh");
    stashFile("$path/haplotypeReads.sh");


  finishStage:
    emitStage($asm, "hap-haplotypeConfigure");

  allDone:
}



sub haplotypeCheck($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $bin     = getBinDirectory();

    my $path    = "haplotype/2-haplotype";

    goto allDone   if (checkHaplotypeReads($asm, "haplotype") == 1);
    goto allDone   if (skipStage($asm, "hap-haplotypeCheck", $attempt) == 1);
    goto allDone   if (sequenceFileExists("$asm.haplotypeReads"));

    #  Compute the size of gkpStore for staging

    setGlobal("corStageSpace", getSizeOfGatekeeperStore($asm));

    #  Figure out if all the tasks finished correctly.

    fetchFile("$path/haplotypeReads.sh");

    my ($jobs, undef) = computeNumberOfHaplotypeJobs($asm);

    my $currentJobID = "0001";
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if (fileExists("$path/results/$currentJobID.success")) {
            push @successJobs, "$path/results/$currentJobID.success\n";

        } else {
            $failureMessage .= "--   job $path/results/$currentJobID.success FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Read haplotype jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Read haplotype jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "hap-haplotypeReads", $attempt);

        submitOrRunParallelJob($asm, "cor", $path, "haplotypeReads", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " read haplotype output files.\n";

    open(L, "> $path/hapjob.files") or caExit("failed to open '$path/hapjob.files'", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/hapjob.files");

    emitStage($asm, "hap-haplotypeReadsCheck");

  allDone:
}




sub dumpHaplotypeReads ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $path    = "haplotype/2-haplotype";
    my $base    = "haplotype";
    my $cmd;

    my @haplotypes = getHaplotypes($base);
    push @haplotypes, "unknown";

    goto allDone   if (checkHaplotypeReads($asm, "haplotype") == 1);
    goto allDone   if (skipStage($asm, "hap-dumpHaplotypeReads") == 1);

    print STDERR "-- Concatenating reads output.\n";

    my $files = 0;
    my $reads = 0;

    stashFile("$path/hapjob.files");

    foreach my $haplotype (@haplotypes) {
       open(O, "| gzip -1c > $asm.${haplotype}Reads.fasta.gz") or caExit("can't open '$asm.${haplotype}Reads.fasta.gz' for writing: $!", undef);
       open(L, "> $asm.${haplotype}Reads.length")              or caExit("can't open '$asm.${haplotype}Reads.length' for writing: $!", undef);

       open(F, "< $path/hapjob.files")                      or caExit("can't open '$path/hapjob.files' for reading: $!", undef);
       open(N, "< $asm.gkpStore/readNames.txt")             or caExit("can't open '$asm.gkpStore/readNames.txt' for reading: $!", undef);


        while (<F>) {
           chomp;

           if (m/^(.*)\/results\/(\d+).success$/) {
              $_ = "$1/results/$2.$haplotype.fasta"
           }
           fetchFile($_);

           open(R, "< $_") or caExit("can't open haplotype output '$_' for reading: $!\n", undef);

           my $h;   #  Current header line
           my $s;   #  Current sequence
           my $n;   #  Next header line

           my $nameid;  #  Currently loaded id and name from gkpStore/readNames
           my $name;

           $n = <R>;  chomp $n;  #  Read the first line, the first header.

           while (!eof(R)) {
              $h = $n;              #  Read name.
               $s = undef;           #  No sequence yet.
               $n = <R>;  chomp $n;  #  Sequence, or the next header.

               #  Read sequence until the next header or we stop reading lines.  Perl seems to be
               #  setting EOF when the last line is read, which is early, IMHO.  We loop until
               #  we both are EOF and have an empty line.

               while (($n !~ m/^>/) && ((length($n) > 0) || !eof(R))) {
                   $s .= $n;

                   $n = <R>;  chomp $n;
               }

               #  Parse the header of the corrected read to find the IID and the split piece number

               my $rid = undef;  #  Read ID

               if ($h =~ m/read(\d+)/) {
                   $rid = $1;
               }

               #  Load the next line from the gatekeeper ID map file until we find the
               #  correct one.

               while (!eof(N) && ($rid != $nameid)) {
                   my $rn = <N>;

                   ($rn =~ m/^(\d+)\s+(.*)$/);

                   $nameid = $1;
                   $name   = $2;
               }

               #  If a match, replace the header with the actual read id.  If no match, use the bogus
               #  corrected read name as is.

               if ($rid eq $nameid) {
                   $h = ">$name id=${rid}";
               }

               #  And write the read to the output as FASTA.

               print O "$h", "\n";
               print O "$s", "\n";

               print L "$name\t", length($s), "\n";

               $reads++;
           }

           $files++;
           close(R);
       }

       close(O);
       close(F);

       stashFile("$asm.${haplotype}Reads.fasta.gz");
       stashFile("$asm.${haplotype}Reads.length");
    }

    #  Now that all outputs are (re)written, cleanup the job outputs.

    print STDERR "--\n";
    print STDERR "-- Purging haplotypeReads output after merging to final output file.\n";
    print STDERR "-- Purging disabled by saveReadHaplotypes=true.\n"  if (getGlobal("saveReadHaplotypes") == 1);

    my $Nsuccess = 0;
    my $Nerr     = 0;
    my $Nfasta   = 0;
    my $Nhap     = 0;

    if (getGlobal("saveReadHaplotypes") != 1) {
        open(F, "< $path/hapjob.files") or caExit("can't open '$path/hapjob.files' for reading: $!", undef);
        while (<F>) {
            chomp;

            if (m/^(.*)\/results\/0*(\d+).success$/) {
                my $ID6 = substr("00000" . $2, -6);
                my $ID4 = substr("000"   . $2, -4);
                my $ID0 = $2;

                if (-e "$1/results/$ID4.success") {
                    $Nsuccess++;
                    unlink "$1/results/$ID4.success";
                }
                if (-e "$1/results/$ID4.err") {
                    $Nerr++;
                    unlink "$1/results/$ID4.err";
                }
                if (-e "$1/results/$ID4.fasta") {
                    $Nfasta++;
                    unlink "$1/results//$ID4.fasta";
                }
                if (-e "$1/results/$ID4.fastaidx") {
                    $Nfasta++;
                    unlink "$1/results/$ID4.fastaidx";
                }

                foreach my $haplotype (@haplotypes) {
                   if (-e "$1/results/$ID4.$haplotype") {
                       $Nhap++;
                       unlink "$1/results/$ID4.$haplotype";
                    }
                    if (-e "$1/results/$ID4.$haplotype.fasta") {
                       $Nhap++;
                       unlink "$1/results/$ID4.$haplotype.fasta";
                    }
                }

            } else {
                caExit("unknown correctReads job name '$_'\n", undef);
            }
        }
        close(F);

        print STDERR "-- Purged $Nsuccess .success sentinels.\n"   if ($Nsuccess > 0);
        print STDERR "-- Purged $Nfasta .fasta inputs.\n"          if ($Nfasta > 0);
        print STDERR "-- Purged $Nerr .err outputs.\n"             if ($Nerr > 0);
        print STDERR "-- Purged $Nhap .out job fasta outputs.\n"   if ($Nhap > 0);
    }

    foreach my $haplotype (@haplotypes) {
        print STDERR "--\n";
        print STDERR "-- Purging haplotype $haplotype gatekeeper store used for classification.\n";
        unlink("$base/$haplotype.gkpStore");
        remove_tree("$haplotype.gkpStore");
     }

  finishStage:
    emitStage($asm, "hap-dumpHaplotypeReads");

  allDone:
    foreach my $haplotype (@haplotypes) {
       print STDERR "--\n";
       print STDERR "-- Haplotype reads saved in '", sequenceFileExists("$asm.${haplotype}Reads"), "'.\n";
    }

    stopAfter("readHaplotyping");
}
