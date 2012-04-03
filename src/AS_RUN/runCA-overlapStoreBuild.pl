#!/usr/bin/env perl

use strict;

use Config;  #  for @signame
use FindBin;
use Cwd;
use Carp;
use FileHandle;

use Sys::Hostname;

use POSIX qw(ceil floor sys_wait_h);

#  usage:  overlapStoreBuild -d $WRK -p $ASM -t obt -g $GKP -i $INP
#
#  Finds overlaps (*.ovb.gz files) in $INP (or any subdirectory), and
#  creates overlap store $WRK/$ASM.${TYP}Store

my $wrk = undef;  #  Path to where the store should be created.
my $asm = undef;  #  Name of our assembly.
my $typ = undef;  #  Type of store to build (obt, dup, mer, ovl).
my $gkp = undef;  #  Path to the gkpStore.
my $inp = undef;  #  Path to input files.

my $bin = getBinDirectory();
my $sbn = getBinDirectoryShellCode();

my $err = 0;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        $wrk = "$ENV{'PWD'}/$wrk" if ($wrk !~ m!^/!);

    } elsif ($arg eq "-p") {
        $asm = shift @ARGV;

    } elsif ($arg eq "-t") {
        $typ = shift @ARGV;

    } elsif ($arg eq "-g") {
        $gkp = shift @ARGV;
        $gkp = "$ENV{'PWD'}/$gkp" if ($gkp !~ m!^/!);

    } elsif ($arg eq "-i") {
        $inp = shift @ARGV;
        $inp = "$ENV{'PWD'}/$inp" if ($inp !~ m!^/!);

    } else {
        die "Unknown option '$arg'\n";
    }
}

$err++  if (!defined($wrk) || (! -d $wrk));
$err++  if (!defined($asm));
$err++  if (!defined($typ));
$err++  if (!defined($gkp) || (! -d $gkp));
$err++  if (!defined($inp) || (! -d $inp));

if ($err) {
    print STDERR "usage: $0 []\n";
    print STDERR "  -d wrk          path to location where store should be created\n";
    print STDERR "  -p asm          prefix of store\n";
    print STDERR "  -t typ          type of store: obt dup mer ovl\n";
    print STDERR "  -g gkp          path to gkpStore\n";
    print STDERR "  -i inp          path to input files\n";
    exit(1);
}






#  Check that all -- or plausibly all -- overlap files are here

my $firstIdx = "999999";
my $lastIdx  = "000000";
my $numJobs  = 0;

my @jobArray;

my $numJobs  = 0;

open(F, "find $inp -type f -name \*.ovb.gz -print | sort |");
while (<F>) {
    chomp;

    if (m/\d\d\d\/(\d\d\d\d\d\d).ovb.gz$/) {
        $firstIdx = $1  if ($1 < $firstIdx);
        $lastIdx = $1   if ($lastIdx < $1);

        push @jobArray, $_;

        $numJobs++;
    }
}
close(F);

$numJobs = scalar(@jobArray);

print STDERR "Found $numJobs jobs from index $firstIdx to $lastIdx.\n";

die "No jobs.\n" if ($numJobs == 0);

#  Submit parallel jobs for bucketizing.  This should really be part
#  of overlap computation itself.

open(F, "> $wrk/1-bucketize.sh") or die;
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
print F "jn=\"undefined\"\n";
print F "\n";

my $tstid = 1;
foreach my $jobname (@jobArray) {
    print F "if [ \"\$jobid\" -eq \"$tstid\" ] ; then jn=\"$jobname\"; fi\n";
    $tstid++;
}

print F "\n";
print F "if [ \$jn = \"undefined\" ] ; then\n";
print F "  echo \"Job out of range.\"\n";
print F "  exit\n";
print F "fi\n";
print F "\n";
print F "if [ -e \"\$jn.bucketizing.gz\" ] ; then\n";
print F "  echo \"Input \$jn is in progress.\"\n";
print F "  exit\n";
print F "fi\n";
print F "\n";
print F "if [ -e \"\$jn.success\" ] ; then\n";
print F "  echo \"Input \$jn finished successfully.\"\n";
print F "  exit\n";
print F "fi\n";
print F "\n";
print F "mv -f \"\$jn\" \"\$jn.bucketizing.gz\"\n";
print F "\n";
print F "$sbn\n";
print F "\n";
print F "\$bin/overlapStoreBucketizer \\\n";
print F "  -o $wrk/$asm.${typ}Store \\\n";
print F "  -g $gkp \\\n";
print F "  -F 512 \\\n";
print F "  -obt \\\n";
print F "  -job \$jobid \\\n";
print F "  -i   \$jn.bucketizing.gz \\\n";
print F "  -e   0.045 \\\n";
print F "&& \\\n";
print F "mv -f \$jn.bucketizing.gz \$jn.success\n";
close(F);


open(F, "> $wrk/2-sort.sh") or die;
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
print F "$sbn\n";
print F "\n";
print F "\$bin/overlapStoreSorter \\\n";
print F "  -o $wrk/$asm.${typ}Store \\\n";
print F "  -g $gkp \\\n";
print F "  -F 512 \\\n";
print F "  -job \$jobid $lastIdx\n";
print F "\n";
print F "\n";
close(F);


open(F, "> $wrk/3-index.sh") or die;
print F "#!/bin/sh\n";
print F "\n";
print F "$sbn\n";
print F "\n";
print F "\$bin/overlapStoreIndexer \\\n";
print F "  -o $wrk/$asm.${typ}Store \\\n";
print F "  -g $gkp \\\n";
print F "\n";
print F "\n";
close(F);


print STDERR "qsub -N ovs1$asm                    -b n -t 1-$numJobs -l memory=1g -j y -o $wrk/1-bucketize-\\\$TASK_ID.err $wrk/1-bucketize.sh\n";
print STDERR "qsub -N ovs2$asm -hold_jid ovs1$asm -b n -t 1-512      -l memory=1g -j y -o $wrk/2-sort-\\\$TASK_ID.err      $wrk/2-sort.sh\n";
print STDERR "qsub -N ovs3$asm -hold_jid ovs2$asm -b n               -l memory=1g -j y -o $wrk/3-index.err                 $wrk/3-index.sh\n";



#
#  Copied from runCA.pl
#

sub getInstallDirectory () {
    my @t = split '/', "$FindBin::RealBin";
    pop @t;                      #  bin
    pop @t;                      #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler

    return($installDir);
}

sub getBinDirectory () {
    my $installDir = getInstallDirectory();

    my $syst = `uname -s`;    chomp $syst;  #  OS implementation
    my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
    my $name = `uname -n`;    chomp $name;  #  Name of the system

    $arch = "amd64"  if ($arch eq "x86_64");
    $arch = "ppc"    if ($arch eq "Power Macintosh");

    return("$installDir/$syst-$arch/bin");
}

sub getBinDirectoryShellCode () {
    my $installDir = getInstallDirectory();
    my $string;

    $string  = "\n";
    $string .= "syst=`uname -s`\n";
    $string .= "arch=`uname -m`\n";
    $string .= "name=`uname -n`\n";
    $string .= "\n";
    $string .= "if [ \"\$arch\" = \"x86_64\" ] ; then\n";
    $string .= "  arch=\"amd64\"\n";
    $string .= "fi\n";
    $string .= "if [ \"\$arch\" = \"Power Macintosh\" ] ; then\n";
    $string .= "  arch=\"ppc\"\n";
    $string .= "fi\n";
    $string .= "\n";
    $string .= "bin=\"$installDir/\$syst-\$arch/bin\"\n";
    $string .= "\n";

    return($string);
}
