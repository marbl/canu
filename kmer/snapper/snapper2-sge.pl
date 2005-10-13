#!/usr/bin/perl

#  Runs snapper2 on SGE, splitting both the genome and query sequences.

use FindBin;

my $genome = "";
my $query  = "";
my $dir    = "";
my $mask   = "";

my $bin = "$FindBin::Bin";

while (scalar(@ARGV) > 0) {
    $arg = shift @ARGV;

    if      ($arg =~ m/^-genome/) {
        $genome = shift @ARGV;
    } elsif ($arg =~ m/^-query/) {
        $query  = shift @ARGV;
    } elsif ($arg =~ m/^-dir/) {
        $dir    = shift @ARGV;
    } elsif ($arg =~ m/^-mask/) {
        $mask   = shift @ARGV;
    } else {
        print STDERR "unknown option '$arg'\n";
    }
}

if (!defined($genome) || !defined($query) || !defined($dir)) {
    print STDERR "usage: $0 <-genome x.fasta> <-query q.fasta> <-dir /path/to/work> <...>\n";
    exit(1);
}

die "Can't find genome '$genome'\n" if (! -e $genome);
die "Can't find queries '$query'\n"  if (! -e $query);

system("mkdir $dir")      if (! -d $dir);
die "Can't find '$dir'\n" if (! -d $dir);

if (! -e "$dir/gen/gen.partitioned") {
    system("mkdir $dir/gen")  if (! -d "$dir/gen");
    system("$bin/leaff -F $genome --partition $dir/gen/gen 2mbp");
    open(F, "> $dir/gen/gen.partitioned");
    close(F);
}

if (! -e "$dir/qry/qry.partitioned") {
    system("mkdir $dir/qry")  if (! -d "$dir/qry");
    system("$bin/leaff -F $query  --partition $dir/qry/qry 32");
    open(F, "> $dir/qry/qry.partitioned");
    close(F);
}

#  Build indexes for everyone -- this prevents the grid jobs from
#  racing to build them.  And it lets us count how many jobs to
#  submit.
#
my $gen = 0;
my $qry = 0;
open(F, "ls $dir/gen/gen-*.fasta $dir/qry/qry-*.fasta |");
while (<F>) {
    chomp;

    if (! -e "${_}idx") {
        system("$bin/leaff -F $_");
    }

    if (m/\/gen-\d\d\d.fasta$/) {
        $gen++;
    } elsif (m/\/qry-\d\d\d.fasta$/) {
        $qry++;
    } else {
        print STDERR "ERROR:  Unknown file '$_'\n";
    }
}
close(F);

if (-e "$mask") {
    system("ln -s $mask $dir/skip.mers.fasta");
}

if (! -e "$dir/skip.mers.fasta") {
    print STDERR "Building skip.mers.fasta!\n";

    system("$bin/meryl -v -B -C -L 100 -m 28 -s $genome -o $dir/gen");
    system("$bin/meryl -Dt -n 200 -s $dir/gen > $dir/skip.mers.fasta");
}

open(F, "> $dir/run.sh");
print F "#!/bin/sh\n";
print F "PIECE=`expr \$SGE_TASK_ID - 1`\n";
print F "GPIECE=`expr \$PIECE % $gen + 1`\n";
print F "QPIECE=`expr \$PIECE / $gen + 1`\n";
print F "GPIECE=`printf %03d \$GPIECE`\n";
print F "QPIECE=`printf %03d \$QPIECE`\n";
print F "echo SGE_TASK_ID = \$SGE_TASK_ID\n";
print F "echo GPIECE      = \$GPIECE\n";
print F "echo QPIECE      = \$QPIECE\n";
#print F "touch $dir/map-gen\$GPIECE-qlt\$QPIECE.success\n";
#print F "exit\n";
print F "$bin/snapper2 -verbose \\\n";
print F "  -mersize 28 -merskip 7 \\\n";
print F "  -genomic $dir/gen/gen-\$GPIECE.fasta \\\n";
print F "  -queries $dir/qry/qry-\$QPIECE.fasta \\\n";
print F "  -mask    $dir/skip.mers.fasta \\\n";
print F "  -output  $dir/map-gen\$GPIECE-qlt\$QPIECE.sim4db \\\n";
print F "  -numthreads 2 \\\n";
print F "  -minmatchidentity 96 \\\n";
print F "  -minmatchcoverage 96 \\\n";
print F "&& \\\n";
print F "touch $dir/map-gen\$GPIECE-qlt\$QPIECE.success\n";
close(F);

my $numJobs = $gen * $qry;

print STDOUT "qsub -t 1-$numJobs -p -50 -j y -o $dir/map-\\\$TASK_ID $dir/run.sh\n";
