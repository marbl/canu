#!/usr/bin/env perl

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
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

my $prefix;
my $FILE     = undef;
my $REF      = undef;
my $withSGE  = 0;
my $withPart = 0;

if (scalar(@ARGV) == 0) {
    die "usage: $0 [-sge] [-partition] PREFIX.tigStore REFERENCE.fasta\n";
}

while (scalar(@ARGV) > 0) {
    if ($ARGV[0] eq "-sge") {
        $withSGE = 1;

    } elsif ($ARGV[0] eq "-partition") {
        $withPart = 1;

    } else {
        $FILE = $ARGV[0];
        $REF  = $ARGV[1];
        shift @ARGV;
    }

    shift @ARGV;
}

die "Must supply 'PREFIX.tigStore' as first non-option argument.\n" if (!defined($FILE));
die "Must supply 'REFERENCE.fasta' as second non-option argument.\n" if (!defined($REF));

if ($FILE =~ m/^(.*).tigStore/) {
    $FILE = $1;
}

if ($FILE =~ m/^(.*)\.\d\d\d\./) {
    $prefix = $1;
} else {
    die "Failed to find prefix in '$FILE', expecting PREFIX.###.tigStore.\n";
}

die "Failed to find reference '$REF'\n" if (! -e $REF);

my %detected1;
my %detected2;


if (! -e "$FILE.tigStore") {
    my $cmd;

    $cmd  = "bogart \\n";
    $cmd .= "  -O ../$prefix.ovlStore \\n";
    $cmd .= "  -G ../$prefix.gkpStore \\n";
    $cmd .= "  -T ../$prefix.tigStore \\n";
    $cmd .= "  -B 75000 \\n";
    $cmd .= "  -e 0.03 \\n";
    $cmd .= "  -E 0 \\n";
    $cmd .= "  -b \\n";
    $cmd .= "  -m 7 \\n";
    $cmd .= "  -U \\n";
    $cmd .= "  -o $prefix \\n";
    $cmd .= "  -D all \\n";
    $cmd .= "  -d overlapQuality \\n";

    print "$cmd\n";

    die "Will not run command.\n";
}

########################################

if ($withPart) {
    my $cmd;

    $cmd .= "gatekeeper";
    $cmd .= " -P $FILE.tigStore.partitioning ../$prefix.gkpStore";

    system($cmd);

    print STDERR "Partitioned.\n";
    exit;
}

########################################

print STDERR "Running consensus\n";
open(R, "> $FILE.sge.sh") or die;
print R "#!/bin/sh\n";
print R "idx=`printf %03d \$SGE_TASK_ID`\n";
print R "utgcns -f -g ../$prefix.gkpStore -t $FILE.tigStore 1 \$idx > $FILE.tigStore/seqDB.v002.p\${idx}.utg.err 2>&1\n";
close(R);

open(R, "> $FILE.tmp.sh") or die;
open(F, "ls $FILE.tigStore/seqDB.v001.p???.dat |");
while (<F>) {
    chomp;

    if (m/seqDB.v001.p(\d\d\d).dat$/) {
        my $idx = $1;

        if (! -e "$FILE.tigStore/seqDB.v002.p${idx}.dat") {
            print STDERR " $_\n";
            print R "utgcns -f -g ../$prefix.gkpStore -t $FILE.tigStore 1 $idx > $FILE.tigStore/seqDB.v002.p${idx}.utg.err 2>&1 &\n";
        } else {
            print STDERR " $_ finished.\n";
        }
    } else {
        die "No '$_'\n";
    }
}
close(F);

print R "wait\n";
close(R);

if ($withSGE) {
    #system("qsub -cwd -j y -o /dev/null -l memory=2g -t 1-9 -q servers.q,black.q $FILE.sge.sh");
    exit;
} else {
    system("sh $FILE.tmp.sh");
}
#unlink "$FILE.tmp.sh";
#unlink "$FILE.sge.sh";

########################################

print STDERR "Dumping fasta\n";
if (! -e "$FILE.fasta") {
    my $cmd;

    $cmd .= "tigStore";
    $cmd .= " -g ../$prefix.gkpStore";
    $cmd .= " -t $FILE.tigStore 2";
    $cmd .= " -U -d consensus";
    $cmd .= "> $FILE.fasta";

    system($cmd);
}

print STDERR "Dumping layout\n";
if (! -e "$FILE.layout") {
    my $cmd;

    $cmd .= "tigStore";
    $cmd .= " -g ../$prefix.gkpStore";
    $cmd .= " -t $FILE.tigStore 2";
    $cmd .= " -U -d layout";
    $cmd .= "> $FILE.layout";

    system($cmd);
}

########################################

print STDERR "Running nucmer\n";
if (! -e "$FILE.delta") {
    my $cmd;

    $cmd .= "nucmer";
    $cmd .= " --maxmatch --coords -p $FILE";
    $cmd .= " $REF";
    $cmd .= " $FILE.fasta";

    system($cmd);
}

if (! -e "$FILE.png") {
    my $cmd;

    $cmd .= "mummerplot";
    $cmd .= " --layout --filter -p $FILE -t png";
    #$cmd .= " --layout          -p $FILE -t png";
    $cmd .= " $FILE.delta";

    system($cmd);
}
