#!/usr/bin/perl

use strict;


my $prefix = "pb1";
my $ctgcns = "../../wgs/Linux-amd64/bin/ctgcns";

my $tv = shift @ARGV;
my $tV = $tv + 1;

die "usage: $0 tv -- recompute tig if it isn't in version 'tv+1'\n"  if (!defined($tv));
die "didn't find '$ctgcns'\n"  if (! -e $ctgcns);

if (! -e "parallel.sh") {
    open(F, "> parallel.sh");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "vv=\$1\n";
    print F "pp=\$2\n";
    print F "id=\$3\n";
    print F "\n";
    print F "export AS_CNS_ERROR_RATE=0.40\n";
    print F "export AS_CGW_ERROR_RATE=0.40\n";
    print F "export AS_OVERLAP_MIN_LEN=250\n";
    print F "\n";
    print F "$ctgcns \\\n";
    print F "  -g ../$prefix.gkpStore -t ../$prefix.tigStore \$vv \$pp \\\n";
    print F "  -c \$id -P 0 -U \\\n";
    print F "  -O parallel.vv\$vv.pt\$pp.id\$id \\\n";
    print F "> parallel.vv\$vv.pt\$pp.id\$id.err 2>&1 \\\n";
    print F "&& \\\n";
    print F "touch parallel.vv\$vv.pt\$pp.id\$id.success\n";
    print F "\n";
    print F "\n";
    close(F);
}

open(CMD, "> submit-parallel-jobs.sh");

for (my $pp=1; $pp < 233; $pp++) {
    my  $vv;

    $tv = sprintf("%03d", $tv);
    $tV = sprintf("%03d", $tV);
    $pp = sprintf("%03d", $pp);

    $vv = $tv if (-e "../$prefix.tigStore/seqDB.v$tv.p$pp.ctg");
    $vv = $tV if (-e "../$prefix.tigStore/seqDB.v$tV.p$pp.ctg");

    last if (!defined($vv));

    print "OPEN version $vv -- ../$prefix.tigStore/seqDB.v$tV.p$pp.ctg\n";

    open(F, "tigStore -g ../$prefix.gkpStore -t ../$prefix.tigStore $vv -cp $pp -D contiglist |");
    while (<F>) {
        my ($maID, $isP, $isD, $ptID, $svID, $fo) = split '\s+', $_;

        next  if ($isD  == 1);
        next  if ($isP  == 0);

        die   if ($svID == 0);

        next  if ($svID > $tv);

        print CMD "qsub -cwd -j y -o /dev/null -b n parallel.sh $tv $pp $maID\n";
    }
    close(F);
}

close(CMD);
