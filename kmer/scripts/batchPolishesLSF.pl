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
#$richUncle = "00006:MRNA:L";
#$richUncle = "00016:BIOINFORMATIC_SVCS:L";

my $numCPU  = 4;
my $sim4db  = "/work/assembly/walenzbp/projects/sim4db/sim4db";
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

    open(F, "> $p.scr");

    for (my $k=0; $k < $size && !eof(IN); $k++) {
        $_ = <IN>;
        print F $_;
    }
    close(F);

    my $cmd;
    $cmd  = "$sim4db -v -YN -aligns -mini 95 -minc 50 ";
    $cmd .= "-cdna    /dev5/walenz/dros/matedFrags.fasta ";
    #$cmd .= "-genomic /dev5/walenz/dros/na_arms.dros.RELEASE3.fasta ";
    $cmd .= "-genomic /dev5/walenz/dros/dros.fasta ";
    $cmd .= "-touch   $p.touch ";
    $cmd .= "-stats   $p.stats ";
    $cmd .= "-output  $p.polished ";
    $cmd .= "-scr     $p.scr ";
    $cmd .= " > $p.answers ";

    print "$cmd\n";

    #if ($numCPU == 0) {
    #    print "No farm!\n";
    #    #system("$bsub $cmd");
    #} else {
    #    &libBri::schedulerSubmit($cmd);
    #    #&libBri::schedulerFinish($numCPU);
    #}
    
    $p++;
}

close(IN);

&libBri::schedulerFinish();
