#!/usr/local/bin/perl

use strict;

if (scalar(@ARGV) < 3) {
    print STDERR "usage: $0 <work-directory> <cdna-file> <genome-file> [sim4 options] < script\n";
    exit;
}

my $sim4db  = "/work/assembly/walenzbp/projects/sim4db/sim4db ";

my $numCPU  = 4;
my $size    = 500000;

my $dir     = shift @ARGV;
my $cdna    = shift @ARGV;
my $genomic = shift @ARGV;

system("mkdir $dir") if (! -d "$dir");

if (! -e "$cdna") {
    print STDERR "Can't find cdna-file '$cdna'\n";
    exit;
}

if (! -e "$genomic") {
    print STDERR "Can't find genome-file '$cdna'\n";
    exit;
}

my $p = "0000";
while (!eof(STDIN)) {

    if (! -e "$dir/$p.sh") {
        print STDERR "Starting $p\n";

        open(F, "| sort -T. -k5n -k3n > $dir/$p.scr");
        for (my $k=0; $k < $size && !eof(STDIN); $k++) {
            $_ = <STDIN>;
            print F $_;
            print Z $_;
        }
        close(F);

        my $cmd = "";
        $cmd .= "$sim4db -YN ";
        $cmd .= "-minidentity 97 -mincoverage 95 -align ";
        $cmd .= "-cdna $cdna ";
        $cmd .= "-genomic $genomic ";
        $cmd .= "-script $dir/$p.scr ";
        $cmd .= "-output $dir/$p.polished ";
        $cmd .= "-stats  $dir/$p.stats ";
        $cmd .= "-touch  $dir/$p.touch ";
        $cmd .= "> $dir/$p.answers";

        open(Z, "> $dir/$p.sh");
        print Z "$cmd\n";
        close(Z);

        system("chmod 755 $dir/$p.sh");

        #system("bsub -q high_mem -o $dir/$p.lsfoutput -R \"select[mem>1300]\" -P 00008:HUM_ANNOT $dir/$p.sh");
    } else {
        print STDERR "Skipping $p\n";

        for (my $k=0; $k < $size && !eof(STDIN); $k++) {
            $_ = <STDIN>;
            print Z $_;
        }
    }

    $p++;
}

#close(STDIN);
