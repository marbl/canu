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

my $jobs        = 512;
my $maxMemory   = 2;    #  gigabytes
my $delete      = "late";

my $maxError    = undef;

my $submit = 0;

my $bin = getBinDirectory();
my $sbn = getBinDirectoryShellCode();

my $err = 0;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg eq "-d") {
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

    } elsif ($arg eq "-jobs") {
        $jobs = shift @ARGV;

    } elsif ($arg eq "-memory") {
        $maxMemory = shift @ARGV;

    } elsif ($arg eq "-deleteearly") {
        $delete = "early";

    } elsif ($arg eq "-nodelete") {
        $delete = "no";

    } elsif ($arg eq "-maxerror") {
        $maxError = shift @ARGV;

    } elsif ($arg eq "-submit") {
        $submit = 1;

    } else {
        die "Unknown option '$arg'\n";
    }
}

if (!defined($wrk) || (! -d $wrk)) {
    print STDERR "ERROR:  Work directory '$wrk' (-d option) not supplied or not found.\n";
    $err++;
}
if (!defined($asm)) {
    print STDERR "ERROR:  Assembly prefix (-p option) not supplied.\n";
    $err++;
}
if (!defined($typ)) {
    print STDERR "ERROR:  Store type (-t option) not supplied.\n";
    $err++;
}
if (($typ ne "obt") &&
    ($typ ne "dup") &&
    ($typ ne "mer") &&
    ($typ ne "ovl")) {
    $err++;
    print STDERR "ERROR:  Store type (-t option) not valid - must be 'obt', 'dup', 'mer' or 'ovl'.\n";
}
if (!defined($gkp) || (! -d $gkp)) {
    print STDERR "ERROR:  Gatekeeper path '$gkp' (-g option) not supplied or not found.\n";
    $err++;
}
if (!defined($inp) || (! -d $inp)) {
    print STDERR "ERROR:  Input path '$inp' (-i option) not supplied or not found.\n";
    $err++;
}

if ($err) {
    print STDERR "\n";
    print STDERR "usage: $0 -d output-directory -p prefix -t type -g gkpStore -i input-directory [options]]\n";
    print STDERR "\n";
    print STDERR "Mandatory parameters:\n";
    print STDERR "  -d wrk          path to location where store should be created\n";
    print STDERR "  -p asm          prefix of store\n";
    print STDERR "  -t typ          type of store: obt dup mer ovl\n";
    print STDERR "  -g gkp          path to gkpStore\n";
    print STDERR "  -i inp          path to input files\n";
    print STDERR "\n";
    print STDERR "Options:\n";
    print STDERR "  -jobs j         create 'j' sorting jobs\n";
    print STDERR "  -memory m       request 'm' gigabytes memory from SGE for sorting\n";
    print STDERR "  -deleteearly    delete intermediate files as soon as possible (unsafe)\n";
    print STDERR "  -nodelete       do not delete intermediate files\n";
    print STDERR "  -maxerror e     discard overlaps with more than 'e' fraction error (no filtering by default)\n";
    print STDERR "  -submit         automatically submit to SGE\n";
    print STDERR "\n";
    print STDERR "This will create an overlap store in three phases.\n";
    print STDERR "\n";
    print STDERR "The first phase will read overlapper outputs found in the '-i inp' directory and write them to the\n";
    print STDERR "store directory.  The store directory is in 'wrk/asm.typStore'.  Overlaps are organized into\n";
    print STDERR "'buckets' with 'slices'.  A bucket corresponds to a single overlapper output, while a slice is a\n";
    print STDERR "range of overlaps (based on fragment ID).  There will be one job for each overlap output file.\n";
    print STDERR "\n";
    print STDERR "The '-jobs j' parameter tells how many slices to make.  More slices means more sort jobs, but each\n";
    print STDERR "job needs less memory.  The maximum number of jobs is dictated by the operating system, as the\n";
    print STDERR "number of open files per process.\n";
    print STDERR "\n";
    print STDERR "The second phase will read a single slice (from multiple buckets) into memory, sort them completely,\n";
    print STDERR "and write the overlaps into the store.  The '-memory m' parameter does NOT control how much memory\n";
    print STDERR "is used, but ONLY tells SGE how much memory each job requires.  There is NO control over how much\n";
    print STDERR "memory is needed, however, if a job needs more memory than requested, it will exit doing nothing.\n";
    print STDERR "\n";
    print STDERR "The third phase will check that all sort jobs in the second phase finished successfully, and create\n";
    print STDERR "a master index for the store.  It is sequential, low memory and generally quick.\n";
    print STDERR "\n";
    print STDERR "Example:\n";
    print STDERR "\n";
    print STDERR "$0 \\\n";
    print STDERR "  -d /work/assembly/godzilla-v1 \\\n";
    print STDERR "  -p godzilla \\\n";
    print STDERR "  -t ovl \\\n";
    print STDERR "  -g /work/assembly/godzilla-v1/godzilla.gkpStore \\\n";
    print STDERR "  -i /scratch/godzillaovl/ \\\n";
    print STDERR "  -jobs 100\n";
    print STDERR "\n";
    print STDERR "This will load overlaps in /scratch/godzillaovl/001, /002, /003, etc into\n";
    print STDERR "store /work/assembly/godzilla-v1/godzilla.ovlStore.\n";
    print STDERR "\n";
    exit(1);
}

#  Check that all -- or plausibly all -- overlap files are here

my $firstIdx = "999999";
my $lastIdx  = "000000";
my $numJobs  = 0;

my @jobArray;

my $numJobs     = 0;

open(F, "find $inp -type f -print |");
while (<F>) {
    chomp;

    if (m/\d\d\d\/(\d\d\d\d\d\d).ovb.[gx]z$/) {
        $firstIdx = $1  if ($1 < $firstIdx);
        $lastIdx = $1   if ($lastIdx < $1);

        push @jobArray, $_;

        $numJobs++;
    }
}
close(F);

@jobArray = sort @jobArray;
$numJobs = scalar(@jobArray);

print STDERR "Found $numJobs jobs from index $firstIdx to $lastIdx.\n";

my $firstNum = int($firstIdx);
my $lastNum  = int($lastIdx);
my $missing  = $lastNum - $firstNum + 1 - $numJobs;

die "First index ($firstIdx) not 000001.\n"  if ($firstIdx ne "000001");

if ($missing != 0) {
    print STDERR "Missing $missing jobs.\n";

    open(F, "> $inp/redo-overlapper.sh") or die;

    my $cc = $firstNum;
    my $ii = 0;

    while ($cc < $lastNum) {
        if ($jobArray[$ii] =~ m/(\d+).ovb/) {
            if ($cc == $1) {
                $ii++;
            } else {
                printf STDERR "  %06d.ovb\n", $cc;
                printf F      "qsub -cwd -j y -o redo-\\\$TASK_ID -pe threads 2 -l mem=2g -t $cc overlap.sh\n";
            }
        } else {
            print "Malformed job '$jobArray[$ii]'\n";
            $ii++;
        }

        $cc++;
    }

    close(F);

    print STDERR "Submission commands written to '$inp/redo-overlapper.sh'\n";

    exit(1);
}

die "No jobs.\n" if ($numJobs == 0);

#  Create an output directory, and populate it with more directories and scripts

system("mkdir -p $wrk/$asm.${typ}Store")                   if (! -d "$wrk/$asm.${typ}Store");
system("mkdir -p $wrk/$asm.${typ}Store/scripts")           if (! -d "$wrk/$asm.${typ}Store/scripts");
system("mkdir -p $wrk/$asm.${typ}Store/logs")              if (! -d "$wrk/$asm.${typ}Store/logs");
system("mkdir -p $wrk/$asm.${typ}Store/logs/1-bucketize")  if (! -d "$wrk/$asm.${typ}Store/logs/bucketize");
system("mkdir -p $wrk/$asm.${typ}Store/logs/2-sort")       if (! -d "$wrk/$asm.${typ}Store/logs/sort");
system("mkdir -p $wrk/$asm.${typ}Store/logs/3-index")      if (! -d "$wrk/$asm.${typ}Store/logs/index");

#  Submit parallel jobs for bucketizing.  This should really be part
#  of overlap computation itself.

open(F, "> $wrk/$asm.${typ}Store/scripts/1-bucketize.sh") or die;
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
print F "bn=`printf %04d \$jobid`\n";
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
print F "if [ -e \"$wrk/$asm.${typ}Store/bucket\$bn/sliceSizes\" ] ; then\n";
print F "  echo \"Bucket $wrk/$asm.${typ}Store/bucket\$bn finished successfully.\"\n";
print F "  exit\n";
print F "fi\n";
print F "\n";
print F "if [ -e \"$wrk/$asm.${typ}Store/create\$bn\" ] ; then\n";
print F "  echo \"Incomplete bucket $wrk/$asm.${typ}Store/create\$bn exists; you must remove to try again.\"\n";
#print F "  rm -rf \"$wrk/$asm.${typ}Store/create\$bn\"\n";
print F "  exit 1\n";
print F "fi\n";
print F "\n";
print F "$sbn\n";
print F "\n";
print F "\$bin/overlapStoreBucketizer \\\n";
print F "  -o $wrk/$asm.${typ}Store \\\n";
print F "  -g $gkp \\\n";
print F "  -F $jobs \\\n";
print F "  -e $maxError \\\n" if (defined($maxError));
print F "  -obt \\\n"   if ($typ eq "obt");
print F "  -dup \\\n"   if ($typ eq "dup");
print F "  -job \$jobid \\\n";
print F "  -i   \$jn\n";
close(F);


open(F, "> $wrk/$asm.${typ}Store/scripts/2-sort.sh") or die;
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
print F "  -delete$delete \\\n"  if ($delete ne "no");
print F "  -M $maxMemory \\\n";
print F "  -o $wrk/$asm.${typ}Store \\\n";
print F "  -g $gkp \\\n";
print F "  -F $jobs \\\n";
print F "  -job \$jobid $lastIdx\n";
print F "\n";
print F "if [ \$? = 0 ] ; then\n";
print F "  echo Success.\n";
print F "else\n";
print F "  echo Failure.\n";
print F "fi\n";
close(F);


open(F, "> $wrk/$asm.${typ}Store/scripts/3-index.sh") or die;
print F "#!/bin/sh\n";
print F "\n";
print F "$sbn\n";
print F "\n";
print F "\$bin/overlapStoreIndexer \\\n";
print F "  -delete \\\n" if ($delete ne "no");
print F "  -o $wrk/$asm.${typ}Store \\\n";
print F "  -F $jobs\n";
print F "\n";
print F "if [ \$? = 0 ] ; then\n";
print F "  echo Success.\n";
print F "else\n";
print F "  echo Failure.\n";
print F "fi\n";
close(F);


my $qsub1 = "";
my $qsub2 = "";
my $qsub3 = "";

$qsub1 .= "qsub -cwd -b n \\\n";
$qsub1 .= "  -N ovs1$asm \\\n";
$qsub1 .= "  -t 1-$numJobs \\\n";
$qsub1 .= "  -l memory=1g \\\n";
$qsub1 .= "  -j y \\\n";
$qsub1 .= "  -o $wrk/$asm.${typ}Store/logs/1-bucketize/bucket-\\\$TASK_ID.err \\\n";
$qsub1 .= "  $wrk/$asm.${typ}Store/scripts/1-bucketize.sh";

$qsub2 .= "qsub -cwd -b n \\\n";
$qsub2 .= "  -N ovs2$asm \\\n";
$qsub2 .= "  -hold_jid ovs1$asm \\\n";
$qsub2 .= "  -t 1-$jobs \\\n";
$qsub2 .= "  -l memory=${maxMemory}g \\\n";
$qsub2 .= "  -j y \\\n";
$qsub2 .= "  -o $wrk/$asm.${typ}Store/logs/2-sort/sort-\\\$TASK_ID.err \\\n";
$qsub2 .= "  $wrk/$asm.${typ}Store/scripts/2-sort.sh";

$qsub3 .= "qsub -cwd -b n \\\n";
$qsub3 .= "  -N ovs3$asm \\\n";
$qsub3 .= "  -hold_jid ovs2$asm \\\n";
$qsub3 .= "  -l memory=1g \\\n";
$qsub3 .= "  -j y \\\n";
$qsub3 .= "  -o $wrk/$asm.${typ}Store/logs/3-index/index.err \\\n";
$qsub3 .= "  $wrk/$asm.${typ}Store/scripts/3-index.sh";

print "$qsub1\n";
system($qsub1) if ($submit == 1);

print "$qsub2\n";
system($qsub2) if ($submit == 1);

print "$qsub3\n";
system($qsub3) if ($submit == 1);


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
