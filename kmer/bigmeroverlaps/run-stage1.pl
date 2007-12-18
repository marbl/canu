#!/usr/bin/perl

use strict;

my $input;
my $dir;

require "runBatchSGE.pl";
require "runCommand.pl";

########################################
sub touch ($) {
    open(F, "> $_[0]");
    close(F);
}
########################################


# perl run-stage1.pl -fasta /home/work/FRAGMENTS/dros/dyakuba.fasta -dir /home/work/bigmertest


while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-fasta/) {
        $input = shift @ARGV;
    } elsif ($arg =~ m/^-dir/) {
        $dir   = shift @ARGV;
    } else {
        die "usage: $0 -fasta some.fasta -dir outputdirectory\n";
    }
}


die "usage: $0 -fasta some.fasta -dir outputdirectory\n" if (!defined($input) || !defined($dir));
die "No input sequences?\n" if (! -e $input);
system("mkdir $dir")      if (! -d $dir);
die "No run directory?\n" if (! -d $dir);

########################################

my $bin;
my $scratch;
my $sge;

$bin     = "/home/work/src/genomics/freebsd/bin";
$bin     = "/bioinfo/assembly/walenz/src/genomics/linux64/bin";

$scratch = "/scratch";

$sge = "-N bigM1";
$sge = "local2";

########################################


sub stage0 {
    print STDERR "stage0\n";

    #
    #  compress the input, build merstream file
    #

    system("mkdir $dir/0-input") if (! -d "$dir/0-input");
    return if (-e "$dir/0-input/frags.success");

    #  XXX  Things that should be global options!
    my $mersize       = 100;

    my $cmd;
    $cmd .= "$bin/markUniqueMers ";
    $cmd .= "  -init ";
    $cmd .= "  -mersize $mersize ";
    $cmd .= "  -f $input ";
    $cmd .= "  -o $dir/0-input/frags";
    if (runCommand($cmd)) {
        rename "$dir/0-input/frags.merStream", 
               "$dir/0-input/frags.merStream.FAILED";
        die "Failed.\n";
    }
    touch("$dir/0-input/frags.success");
}

sub stage1 ($) {
    print STDERR "stage1\n";

    my $pass = shift @_;  #  1 or 2
    my $name = "";

    if ($pass == 1) {
        $name = "1-markDistinct";
    } elsif ($pass == 2) {
        $name = "2-markDistinct";
    } else {
        die "Invalid pass '$pass'\n";
    }

    system("mkdir $dir/$name") if (! -d "$dir/$name");
    return if (-e "$dir/$name/marked.success");

    #  We need to analyze the number of fragments we have, and pick
    #  hashsize appropriately.  Then we can pick partitionbits to fit into
    #  our desired memory space.

    #  Hash sizes:
    #    34 = 16 billion entries = 2048MB data
    #    33 =  8 billion entries = 1024MB data
    #    32 =  4 billion entries =  512MB data

    #  XXX  Things that should be global options!
    my $mersize       = 100;
    my $hashsize      = 34;
    my $partitionbits = 31;
    my $iters         = 2;

    my $seed          = time();

    if (! -e "$dir/$name/mark.sh") {
        open(F, "> $dir/$name/mark.sh");
        print F "#!/bin/sh\n";
        print F "\n";
        print F "jobid=\$SGE_TASK_ID\n";
        print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
        print F "  jobid=\$1\n";
        print F "fi\n";
        print F "if [ x\$jobid = x ]; then\n";
        print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "jobid=`expr \$jobid - 1`\n";
        print F "jobid=`printf %03d \$jobid`\n";
        print F "\n";
        print F "rm -f $scratch/pass$pass-\$jobid.markbits\n";
        print F "\n";
        print F "$bin/markUniqueMers \\\n";
        print F "  -pass1 \\\n"                                         if ($pass == 1);
        print F "  -pass2 ../1-markDistinct/pass1-merged.markbits \\\n" if ($pass == 2);
        print F "  -mersize $mersize \\\n";
        print F "  -hashsize $hashsize \\\n";
        print F "  -seed $seed \\\n";
        print F "  -partitionbits $partitionbits \\\n";
        print F "  -partition \$jobid \\\n";
        print F "  -iters $iters \\\n";
        #print F "  -f $input \\\n";
        print F "  -m $dir/0-input/frags \\\n";
        print F "  -o $scratch/pass$pass-\$jobid.markbits \\\n";
        print F "&& \\\n";
        print F "bzip2 -9 $scratch/pass$pass-\$jobid.markbits \\\n";
        print F "&& \\\n";
        print F "cp -p $scratch/pass$pass-\$jobid.markbits.bz2 \\\n";
        print F "      $dir/$name/pass$pass-\$jobid.markbits.bz2 \\\n";
        print F "&& \\\n";
        print F "touch $dir/$name/pass$pass-\$jobid.markbits.success\n";
        print F "\n";
        print F "rm -f $scratch/pass$pass-\$jobid.markbits.bz2\n";
        print F "rm -f $scratch/pass$pass-\$jobid.markbits\n";
        close(F);

        system("chmod 775 $dir/$name/mark.sh");
    }

    #  Perl doesn't seem to have a pow() function??
    #
    my $jobs = 1;
    for (my $x=0; $x<($hashsize - $partitionbits); $x++) {
        $jobs *= 2;
    }

    if (0 < runBatchSGE("$dir/$name",
                        $jobs,
                        "$dir/$name/pass$pass-###.markbits.success",
                        "$dir/$name/mark.sh",
                        $sge)) {
        print STDERR "Jobs submitted!\n";
        exit(0);
    }

    #  Otherwise, all jobs are done.

    touch("$dir/$name/marked.success");
}


sub stage1merge ($) {
    print STDERR "stage1merge\n";

    my $pass = shift @_;  #  1 or 2
    my $name = "";

    if ($pass == 1) {
        $name = "1-markDistinct";
    } elsif ($pass == 2) {
        $name = "2-markDistinct";
    } else {
        die "Invalid pass '$pass'\n";
    }

    return if (-e "$dir/$name/merged.success");
    die    if (! -e "$dir/$name/marked.success");

    my $cmd;
    $cmd  = "$bin/mergeMarkings \\\n";
    $cmd .= " -o $dir/$name/pass$pass-merged.markbits \\\n";

    open(F, "ls $dir/$name/pass$pass-???.markbits.bz2 |");
    while (<F>) {
        chomp;
        $cmd .= " $_";
    }
    close(F);

    if (runCommand($cmd)) {
        rename "$dir/$name/pass$pass-merged.markbits",
               "$dir/$name/pass$pass-merged.markbits.FAILED";
        die "Failed.\n";
    }

    touch("$dir/$name/merged.success");
}



sub stage3 () {
    print STDERR "stage3\n";

    my $name = "3-propagate";

    system("mkdir $dir/$name") if (! -d "$dir/$name");
    return if (-e "$dir/$name/marked.success");

    #  dros has 2^31 mers.

    my $mersize       = 100;
    my $hashsize      = 38;
    my $partitionbits = 31;

    my $seed          = time();

    if (! -e "$dir/$name/mark.sh") {
        open(F, "> $dir/$name/mark.sh");
        print F "#!/bin/sh\n";
        print F "\n";
        print F "jobid=\$SGE_TASK_ID\n";
        print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
        print F "  jobid=\$1\n";
        print F "fi\n";
        print F "if [ x\$jobid = x ]; then\n";
        print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "jobid=`expr \$jobid - 1`\n";
        print F "jobid=`printf %03d \$jobid`\n";
        print F "\n";
        print F "rm -f $scratch/propagated-\$jobid.markbits\n";
        print F "\n";
        print F "$bin/markUniqueMers \\\n";
        print F "  -propagate ../2-markDistinct/pass2-merged.markbits \\\n";
        print F "  -mersize $mersize \\\n";
        print F "  -hashsize $hashsize \\\n";
        print F "  -seed $seed \\\n";
        print F "  -partitionbits $partitionbits \\\n";
        print F "  -partition \$jobid \\\n";
        print F "  -m $dir/0-input/frags \\\n";
        print F "  -o $scratch/propagated-\$jobid.markbits \\\n";
        print F "&& \\\n";
        print F "bzip2 -9 $scratch/propagated-\$jobid.markbits \\\n";
        print F "&& \\\n";
        print F "cp -p $scratch/propagated-\$jobid.markbits.bz2 \\\n";
        print F "      $dir/$name/propagated-\$jobid.markbits.bz2 \\\n";
        print F "&& \\\n";
        print F "touch $dir/$name/propagated-\$jobid.markbits.success\n";
        print F "\n";
        print F "rm -f $scratch/propagated-\$jobid.markbits.bz2\n";
        print F "rm -f $scratch/propagated-\$jobid.markbits\n";
        close(F);

        system("chmod 775 $dir/$name/mark.sh");
    }

    #  Perl doesn't seem to have a pow() function??
    #
    my $jobs = 1;
    for (my $x=0; $x<($hashsize - $partitionbits); $x++) {
        $jobs *= 2;
    }

    if (0 == runBatchSGE("$dir/$name",
                         $jobs,
                         "$dir/$name/propagated-###.markbits.success",
                         "$dir/$name/mark.sh",
                         $sge)) {
        touch("$dir/$name/marked.success");
    } else {
        print STDERR "Jobs submitted!\n";
        exit(0);
    }
}




stage0();
stage1(1);
stage1merge(1);
stage1(2);
stage1merge(2);
stage3();
