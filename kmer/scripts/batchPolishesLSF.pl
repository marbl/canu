#!/usr/local/bin/perl

$| = 1;

use strict;
use Fcntl;

#  We use the libBri package of usefull stuff.  It's located in the same place
#  as the ESTmapper.pl script.  That's what FindBin tells us.
#
use FindBin;
use lib "$FindBin::Bin";
use libBri;

my $richUncle;
$richUncle = "00006:MRNA:L";
$richUncle = "00016:BIOINFORMATIC_SVCS:L";

my $numCPU  = 4;
my $sim4db  = "/home/walenzbp/projects/ESTmapper/sim4dbseq ";
my $cdna    = "/dev5/walenz/TDW-20020226-dbEST-human/0-input/cDNA.fasta";
my $genomic = "/dev5/walenz/TDW-20020226-dbEST-human/0-input/genomic.fasta";
#my $bsub    = "bsub -q assembly -o /dev/null -R \"select[physmem>400]rusage[physmem=800]\" -P $richUncle ";
my $size    = 100000;

if ($ARGV[0] eq "-lsf") {
    shift @ARGV;
    $numCPU = 0;
}

if ($ARGV[0] eq "-ncpu") {
    shift @ARGV;
    $numCPU = shift @ARGV;

    &libBri::schedulerSetNumberOfProcesses($numCPU);
    &libBri::schedulerSetNumberOfProcessesToWaitFor($numCPU);
} else {
    &libBri::schedulerSetNumberOfProcesses(4);
    &libBri::schedulerSetNumberOfProcessesToWaitFor(0);
}

open(IN, "< $ARGV[0]");

my $p = "0000";
while (!eof(IN)) {

    #  Try top open the script file for exclusive writing.  Only
    #  one process should succeed this.
    #
    if (sysopen(F, "$p.scr", O_WRONLY | O_EXCL | O_CREAT)) {
        print STDERR "Starting $p\n";

        open(F, "> $p.scr");

        for (my $k=0; $k < $size && !eof(IN); $k++) {
            $_ = <IN>;
            print F $_;
        }
        close(F);

        my $cmd;
        $cmd  = "$sim4db -v -cdna $cdna -genomic $genomic ";
        $cmd .= "-script $p.scr ";
        $cmd .= "-output $p.polished ";
        $cmd .= "-stats  $p.stats ";
        $cmd .= "-mincoverage 45 -minidentity 90 -cut 0.6 -align ";

        #print "$cmd\n";

        if ($numCPU == 0) {
            print "No farm!\n";
            #system("$bsub $cmd");
        } else {
            &libBri::schedulerSubmit($cmd);
            #&libBri::schedulerFinish($numCPU);
        }
    } else {
        print STDERR "Skipping $p\n";
            
        for (my $k=0; $k < $size && !eof(IN); $k++) {
            $_ = <IN>;
        }
    }
    
    $p++;
}

close(IN);

&libBri::schedulerFinish();
