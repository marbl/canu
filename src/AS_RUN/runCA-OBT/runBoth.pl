#!/usr/bin/env perl

use FindBin;

my $argS = join ' ', @ARGV;
$argS =~ s/-d (\S+) //;
my $outDir = $1;

my $runCAOut = "${outDir}_runCA";
my $carunOut = "${outDir}_carun";
my $runCA = "$FindBin::Bin/runCA -d $runCAOut";
my $carun = "$FindBin::Bin/carun -d $carunOut";

my $jobId = `$carun $argS | awk '/^Your job / { print \$3 }'`;
chomp $jobId;
print "jobId:$jobId\n";
my $qstat = `qstat -j $jobId | wc -l`;
chomp $qstat;
if ( $qstat < 3 ) {
    die "qstat shows nothing";
}
system "$runCA $argS > runCA.out 2>&1";
while ($qstat > 2 ) {
    sleep 120;
    my $jobs = `grep 'Your job ' $carunOut/runCA.sge.out.* | awk '{ print \$3 }'`;
    $jobs =~ tr/\n/,/;
    $jobs .= $jobId;
    $qstat = `qstat -j $jobs 2>&1 | wc -l`;
    chomp $qstat;
}
print "Comparing qc files\n";
system "cmp $runCAOut/*.qc $carunOut/*.qc";
