#!/usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin/util";
#use scheduler;
#use fasta;

require "run.pl";


#  configureESTmapper.pl
#
#  Builds tables, prepares the genome for ESTmapper.
#
#  What are we doing?
#    -genome    some.fasta
#
#  Where?
#    -genomedir some-directory
#
#  How will we run?  Exactly one of -memory and -segments should be specified.
#    -mersize   k
#    -merskip   s
#    -memory    max memory per segment
#    -segments  number of segments
#
#  Where will we compute what we need to compute?
#    -sge       sge-options, e.g., "-pe thread 2", accounting info, etc.
#    -local

my $exechome  = "$FindBin::Bin";
my $leaff     = "$exechome/leaff";
my $posdb     = "$exechome/positionDB";
my $meryl     = "$exechome/meryl";
my $mimsf     = "$exechome/mersInMerStreamFile";

my $genome    = undef;
my $genomedir = undef;
my $mersize   = 20;
my $merskip   = 0;
my $memory    = 1000;
my $segments  = 0;
my $local     = 1;
my $sge       = undef;
my $sgename   = "EMconfig";

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg eq "-genome") {
        $genome = shift @ARGV;
    } elsif (($arg eq "-genomedir") || ($arg eq "-path")) {
        $genomedir = shift @ARGV;
    } elsif ($arg eq "-mersize") {
        $mersize  = int(shift @ARGV);
    } elsif ($arg eq "-merskip") {
        $merskip  = int(shift @ARGV);
    } elsif ($arg eq "-memory") {
        $memory = int(shift @ARGV);
    } elsif ($arg eq "-segments") {
        $segments = int(shift @ARGV);
    } elsif ($arg eq "-sge") {
        $local = undef;
        $sge   = shift @ARGV;
    } elsif ($arg eq "-sgename") {
        $sgename = shift @ARGV;
    } elsif ($arg eq "-local") {
        $local = 1;
    } elsif ($arg eq "-h") {
        undef $genome;
        undef $genomedir;
        undef @ARGV;
    } elsif ($arg eq "-justtestingifitworks") {
        exit(0);
    } else {
        die "ERROR: unknown arg '$arg'\n";
    }
}
if (!defined($genome) || !defined($genomedir)) {
    print STDERR "usage: $0 -genome g.fasta -path /some/path [args]\n";
    print STDERR "  -genome g.fasta   the genome to map to\n";
    print STDERR "  -genomedir d      the directory to save the configuration in\n";
    print STDERR "\n";
    print STDERR "  -mersize m        use m-mers (default 20)\n";
    print STDERR "  -merskip s        skip s m-mers between mers (default 0, use all mers)\n";
    print STDERR "  -memory M         use M MB memory for the search processes (default 1000MB)\n";
    print STDERR "  -segments S       use S search processes (default, based on memory size)\n";
    print STDERR "  -sge              compute the configuration on the grid; args are passed to qsub\n";
    print STDERR "  -sgename          sge job name (default 'EMconfig')\n";
    print STDERR "  -local            compute the configuration right now (the default)\n";
    print STDERR "\n";
    print STDERR "Example:\n";
    print STDERR "  configureESTmapper.pl -genome B35LC.fasta -genomedir B35LC -memory 900 -sge \"-pe thread 2\"\n";
    print STDERR "\n";
    exit(1);
}

$genome    = "$ENV{'PWD'}/$genome"     if ($genome    !~ m!^/!);
$genomedir = "$ENV{'PWD'}/$genomedir"  if ($genomedir !~ m!^/!);

system("mkdir -p $genomedir") if (! -d $genomedir);

if ($genome !~ m/^\//) {
    my $cwd = `pwd`;
    chomp $cwd;
    $genome = "$cwd/$genome";
}

die "Can't find genome '$genome'\n"         if (! -e $genome);
die "Can't find output directory '$genomedir'\n" if (! -d $genomedir);

print STDERR "Configuring ESTmapper:\n";
print STDERR "  merSize $mersize\n";
print STDERR "  merSkip $merskip\n";
print STDERR "  ${memory}MB\n" if (defined($memory));
print STDERR "  $segments segments\n" if (defined($segments));

symlink "${genome}",    "$genomedir/genome.fasta"    if ((! -f "$genomedir/genome.fasta"));
symlink "${genome}idx", "$genomedir/genome.fastaidx" if ((! -f "$genomedir/genome.fastaidx") && (-f "${genome}idx"));

if (! -f "$genomedir/genome.fastaidx") {
    print STDERR "configureESTmapper-- Generating the genome index.\n";
    if (runCommand("$leaff -F $genomedir/genome.fasta")) {
        unlink "$genomedir/genome.fastaidx";
        die "Failed.\n";
    }
}


#  Build a seqStore so that we can run the positionDB builders in
#  parallel...from here, we can guesstimate an optimal partitioning,
#  though that will be impossible without actually contstructing every
#  mer (consider spaced and/or compressed mers) and so we live with
#  it.
#

print STDERR "configureESTmapper-- Initializing positionDB creation.\n";

if ((! -e "$genomedir/genome.seqStore.blocks") ||
    (! -e "$genomedir/genome.seqStore.sequence") ||
    (! -e "$genomedir/genome.seqStore.out")) {
    if (runCommand("$mimsf $genomedir/genome.fasta $genomedir/genome > $genomedir/genome.seqStore.out 2>&1")) {
        die "Failed.\n";
    }
}

my $acgtInFile     = 0;
my $acgtPerSegment = 0;
my $segmentOverlap = 10000000;

open(F, "< $genomedir/genome.seqStore.out") or die;
while (<F>) {
    if (m/^Found\s+(\d+)\s+ACGT/) {
        $acgtInFile = $1;
    }
}
close(F);

print STDERR "Found $acgtInFile ACGT in the input.\n";
die "No mers found?\n" if ($acgtInFile <= 0);

#  XXX:  Magic Number!  12 bytes per base!

if ($memory > 0) {
    $acgtPerSegment = int($memory / 12 * 1000000) + 1;
    print STDERR "configureESTmapper-- packing to preserve ${memory}MB memory limit ($acgtPerSegment mers per segment)\n";
}

if ($segments > 0) {
    $acgtPerSegment = int($acgtInFile / $segments + $segmentOverlap) + 1;
    print STDERR "configureESTmapper-- packing to preserve $segments processor limit ($acgtPerSegment mers per segment)\n";
}

$memory = int($acgtPerSegment * 12 / 1000000);

open(F, "> $genomedir/memoryLimit") or die "Can't write $genomedir/memoryLimit\n";
print F "$memory\n";
close(F);


my $merBeg = 0;
my $merEnd = 0;
my $segId  = "000";

open(F, "> $genomedir/segments");
open(S, "> $genomedir/create.dat");
while ($merBeg < $acgtInFile) {
    $merEnd  = $merBeg + $acgtPerSegment;

    print F "$segId\n";
    print S "$segId $merBeg $merEnd\n";

    $merBeg += $acgtPerSegment - $segmentOverlap;
    $segId++;
}
close(F);
close(S);

print STDERR "configureESTmapper-- Created $segId groups with maximum memory requirement of ${memory}MB.\n";
die "Created no groups?\n" if (int($segId) == 0);


#  Configure meryl
#
if (! -e "$genomedir/genome.merylArgs") {
    my $cmd;
    $cmd  = "$meryl";
    $cmd .= " -B -L 5 -f -m $mersize -segments $segId -configbatch";
    $cmd .= " -s $genomedir/genome.fasta";
    $cmd .= " -o $genomedir/genome";
    $cmd .= " > $genomedir/meryl.config.out 2>&1";
    if (runCommand($cmd)) {
        die "Failed.\n";
    }
}


#  Create the script that builds the positionDB's and meryl partitions
#
open(F, "> $genomedir/create.sh");
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
print F "jobp=`cat $genomedir/create.dat | head -\$jobid | tail -1`\n";
print F "\n";
print F "seg=`echo \$jobp | awk '{ print \$1 }'`\n";
print F "beg=`echo \$jobp | awk '{ print \$2 }'`\n";
print F "end=`echo \$jobp | awk '{ print \$3 }'`\n";
print F "\n";
print F "if [ ! -e \"$genomedir/genome.seqStore.sequence\" ] ; then\n";
print F "  echo \"Didn't find the seqStore!\"\n";
print F "  exit 1\n";
print F "fi\n";
print F "\n";
print F "ln -s \"$genomedir/genome.seqStore.blocks\"   \"$genomedir/seg\$seg.building.posDB.seqStore.blocks\"\n";
print F "ln -s \"$genomedir/genome.seqStore.sequence\" \"$genomedir/seg\$seg.building.posDB.seqStore.sequence\"\n";
print F "\n";
print F "\n";
print F "if [ ! -e \"$genomedir/seg\$seg.posDB\" ] ; then\n";
print F "  $posdb \\\n";
print F "    -mersize $mersize \\\n";
print F "    -merbegin \$beg \\\n";
print F "    -merend \$end \\\n";
print F "    -sequence \"$genomedir/genome.fasta\" \\\n";
print F "    -output   \"$genomedir/seg\$seg.building.posDB\" \\\n";
print F "  && \\\n";
print F "  mv \"$genomedir/seg\$seg.building.posDB\" \\\n";
print F "     \"$genomedir/seg\$seg.posDB\"\n";
print F "fi\n";
print F "\n";
print F "\n";
print F "$meryl \\\n";
print F "  -countbatch `expr \$jobid - 1` \\\n";
print F "  -v \\\n";
print F "  -o \"$genomedir/genome\" \\\n";
print F "|| \\\n";
print F "rm -f \"$genomedir/genome.mcidx\" \"$genomedir/genome.mcdat\"\n";
print F "\n";
print F "\n";
print F "rm -f \"$genomedir/seg\$seg.building.posDB.seqStore.blocks\"\n";
print F "rm -f \"$genomedir/seg\$seg.building.posDB.seqStore.sequence\"\n";
close(F);


#  Create the script that merges meryl outputs
#
open(F, "> $genomedir/meryl.sh");
print F "#!/bin/sh\n";
print F "\n";
print F "if [ ! -e \"$genomedir/genome.mcidx\" ] ; then\n";
print F "  $meryl -mergebatch -o \"$genomedir/genome\"\n";
print F "fi\n";
print F "if [ ! -e \"$genomedir/frequentMers-ge1000.fasta\" ] ; then\n";
print F "  $meryl -Dt -n 1000 -s \"$genomedir/genome\" \\\n";
print F "    > \"$genomedir/frequentMers-ge1000.fasta\"\n";
print F "fi\n";
close(F);



########################################
#
#  run the jobs.
#
if      ($local) {
    my $seg = "000";

    while ($seg ne $segId) {
        #  Copy $seg (a string) into $s (an integer).
        my $s = int($seg);

        print STDERR "Creating $seg out of $segId\n";
        
        if ((! -e "$genomedir/seg$seg.posDB") || (! -e "$genomedir/genome.batch$s.mcdat")) {
            $s++;
            runCommand("/bin/sh $genomedir/create.sh $s") and die "Segment $seg failed.\n";
        }

        $seg++;
        $seg = substr("000$seg", -3);
    }
    runCommand("/bin/sh $genomedir/meryl.sh") and die "Meryl failed.\n";
} elsif ($sge) {

    #  Check if we need to submit pieces of the array, or if we can submit the whole thing.
    #
    my @ap;
    my $wholeThing = 0;

    system("mkdir $genomedir/sgeout") if (! -d "$genomedir/sgeout");

    my $sgebuildname = "$sgename." . time();

    my $seg = "000";
    while ($seg ne $segId) {
        if (-e "$genomedir/seg$seg.posDB") {
            #print STDERR "Segment $seg finished successfully!\n";
        } else {
            #print STDERR "Segment $seg failed.\n";
            $ap[$seg] = 1;
            $wholeThing++;
        }

        $seg++;
        $seg = substr("000$seg", -3);
    }

    if ($wholeThing == $seg) {
        #  Yippee!  Submit all at once!
        #
        if (runCommand("qsub -cwd -j y -o $genomedir/sgeout/seg\\\$TASK_ID.out -t 1-$segId $sge -N $sgebuildname $genomedir/create.sh")) {
            die "SGE submission failed?\n";
        }
    } elsif ($wholeThing > 0) {
        #  Dang, we need to submit individually....or we can take five
        #  minutes and figure out ranges to submit.
        #
        my $st;
        my $ed;
        my $it = 0;

        #  +2 so that we run off the end -- ensuring that we submit
        #  even the last batch of jobs.

        while ($it < $segId + 2) {
            if (!defined($st) && ($ap[$it] == 1)) {
                #  SGE wants to start at 1, we start at 0.
                $st = $it + 1;
            }
            if (defined($st) && !defined($ed) && ($ap[$it] == 0)) {
                #  SGE wants to start at 1, we start at 0.
                $ed = $it;
            }
            if (defined($st) && defined($ed)) {
                #print STDERR "submit $st - $ed\n";
                if (runCommand("qsub -cwd -j y -o $genomedir/sgeout/seg\\\$TASK_ID.out -t $st-$ed $sge -N $sgebuildname $genomedir/create.sh")) {
                    die "SGE submission failed?\n";
                }
                undef $st;
                undef $ed;
            }

            $it++;
        }
    } else {
        print STDERR "All segments computed successfully!\n";
    }
    if (runCommand("qsub -cwd -j y -o $genomedir/sgeout/meryl.out $sge -hold_jid $sgebuildname -N $sgename $genomedir/meryl.sh")) {
        die "SGE submission failed?\n";
    }
} else {
    die "HELP!  I don't know how to run jobs!\n";
}

