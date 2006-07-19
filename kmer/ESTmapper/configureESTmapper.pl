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
#    -path      some-directory
#
#  How will we run?  Exactly one of -memory and -segments should be specified.
#    -mersize   k
#    -merskip   s
#    -memory    max memory per segment
#    -segments  number of segments
#
#  To get the list of frequent kmers.  If not specified, it will be computed.
#    -species   human | mouse | rat | none
#
#  Where will we compute what we need to compute?
#    -sge       sge-options, e.g., "-pe thread 2", accounting info, etc.
#    -local

my $exechome  = "$FindBin::Bin";
my $leaff     = "$exechome/leaff";
my $posdb     = "$exechome/positionDB";
my $mimsf     = "$exechome/mersInMerStreamFile";

my $genome    = undef;
my $path      = undef;
my $mersize   = 20;
my $merskip   = 0;
my $memory    = 1000;
my $segments  = 0;
my $species   = undef;
my $local     = 1;
my $sge       = undef;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg eq "-genome") {
        $genome = shift @ARGV;
    } elsif ($arg eq "-path") {
        $path = shift @ARGV;
    } elsif ($arg eq "-mersize") {
        $mersize  = int(shift @ARGV);
    } elsif ($arg eq "-merskip") {
        $merskip  = int(shift @ARGV);
    } elsif ($arg eq "-memory") {
        $memory = int(shift @ARGV);
    } elsif ($arg eq "-segments") {
        $segments = int(shift @ARGV);
    } elsif ($arg eq "-species") {
        $species = shift @ARGV;
    } elsif ($arg eq "-sge") {
        $local = undef;
        $sge   = shift @ARGV;
    } elsif ($arg eq "-local") {
        $local = 1;
    } else {
        die "ERROR: unknown arg '$arg'\n";
    }
}

if (!defined($genome) || !defined($path)) {
    print STDERR "usage: $0 -genome g.fasta -path /some/path [args]\n";
    print STDERR "  -genome g.fasta          The genome to map to.\n";
    print STDERR "  -path d                  Absolute path to the directory to save the configuration.\n";
    print STDERR "\n";
    print STDERR "Example:\n";
    print STDERR "  time perl /home/work/src/genomics/ESTmapper/configureESTmapper.pl -genome B35LC.fasta -path /n6/junk -memory 900 -sge \"-pe thread 2\"\n";
    print STDERR "\n";
    exit(1);
}

system("mkdir -p $path") if (! -d $path);

if ($genome !~ m/^\//) {
    my $cwd = `pwd`;
    chomp $cwd;
    $genome = "$cwd/$genome";
}

die "Can't find genome '$genome'\n"         if (! -e $genome);
die "Can't find output directory '$path'\n" if (! -d $path);

print STDERR "Configuring ESTmapper:\n";
print STDERR "  merSize $mersize\n";
print STDERR "  merSkip $merskip\n";
print STDERR "  ${memory}MB\n" if (defined($memory));
print STDERR "  $segments segments\n" if (defined($segments));
print STDERR "  $species\n" if (defined($species));

symlink "${genome}",    "$path/genome.fasta"    if ((! -f "$path/genome.fasta"));
symlink "${genome}idx", "$path/genome.fastaidx" if ((! -f "$path/genome.fastaidx") && (-f "${genome}idx"));

if (! -f "$path/genome.fastaidx") {
    print STDERR "configureESTmapper-- Generating the genome index.\n";
    if (runCommand("$leaff -F $path/genome.fasta")) {
        unlink "$path/genome.fastaidx";
        die "Failed.\n";
    }
}


#  Build the merStreamFile so that we can run the positionDB builders
#  in parallel...it also lets us know exactly how many mers are in the
#  input so we can precisely partition the genome.

print STDERR "configureESTmapper-- Initializing positionDB creation.\n";

if (! -e "$path/init.merStream") {
    my $cmd;
    $cmd  = "$posdb";
    $cmd .= " -mersize $mersize";
    $cmd .= " -merskip $merskip";
    $cmd .= " -merbegin 1 -merend 100";
    $cmd .= " -sequence $path/genome.fasta";
    $cmd .= " -output $path/init";
    if (runCommand($cmd)) {
        die "Failed.\n";
    }

    unlink("$path/init");
}


#  Figure out how many mers are in our sequence, then decide on a
#  partitioning.

#  mersPerSegment - the _actual_ number of mers we can stuff into a
#  segment.  We then overlap each segment by 1,000,000 mers.

my $mersInFile     = `$mimsf $path/init`;  chomp $mersInFile;
my $mersPerSegment = 0;
my $segmentOverlap = 10000000;

#  XXX:  Magic Number!  12 bytes per base!

if ($memory > 0) {
    print STDERR "configureESTmapper-- packing to preserve ${memory}MB memory limit\n";
    $mersPerSegment = int($memory / 12 * 1000000) + 1;
}

if ($segments > 0) {
    print STDERR "configureESTmapper-- packing to preserve $segments processor limit\n";
    $mersPerSegment = int($mersInFile / $segments + $segmentOverlap) + 1;
}

$memory = int($mersPerSegment * 12 / 1000000);

open(F, "> $path/memoryLimit") or die "Can't write $path/memoryLimit\n";
print F "$memory\n";
close(F);



#  Create the script that builds the positionDB's
#
open(F, "> $path/create.sh");
print F "#!/bin/sh\n";
print F "/usr/bin/perl $path/create.pl \$@\n";
close(F);

open(F, "> $path/create.pl");
print F "#!/usr/bin/perl\n";
print F "use strict;\n";
print F "open(F, \"< $path/create.dat\") or die \"Failed to open '$path/create.dat'\\n\";\n";
print F "my \@args = <F>;\n";
print F "close(F);\n";
print F "chomp \@args;\n";
print F "\n";
print F "my \$jid = shift @ARGV;\n";
print F "if (!defined(\$jid)) {\n";
print F "    \$jid = \$ENV{'SGE_TASK_ID'} - 1;\n";
print F "}\n";
print F "\n";
print F "my (\$seg, \$beg, \$end) = split '\\s+', \$args[\$jid];\n";
print F "\n";
print F "if (! -e \"$path/init.merStream\") {\n";
print F "  die \"Didn't find the merStreamFile!\\n\";\n";
print F "}\n";
print F "\n";
print F "system(\"ln -s $path/init.merStream $path/seg\$seg.building.posDB.merStream\");\n";
print F "\n";
print F "my \$cmd;\n";
print F "\$cmd  = \"$posdb -mersize $mersize -merbegin \$beg -merend \$end -sequence $path/genome.fasta -output $path/seg\$seg.building.posDB\";\n";
print F "\$cmd .= \" && mv $path/seg\$seg.building.posDB $path/seg\$seg.posDB\";\n";
print F "print STDERR \"\$cmd\\n\";\n";
print F "\n";
print F "system(\$cmd);\n";
print F "\n";
print F "unlink(\"$path/seg\$seg.building.posDB.merStream\");\n";
close(F);


my $merBeg = 0;
my $merEnd = 0;
my $segId  = "000";

open(F, "> $path/segments");
open(S, "> $path/create.dat");
while ($merBeg < $mersInFile) {
    $merEnd  = $merBeg + $mersPerSegment;

    print F "$segId\n";
    print S "$segId $merBeg $merEnd\n";

    $merBeg += $mersPerSegment - $segmentOverlap;
    $segId++;
}
close(F);
close(S);

print STDERR "configureESTmapper-- Created $segId groups with maximum memory requirement of ${memory}MB.\n";

########################################
#
#  run the jobs.
#
if      ($local) {
    my $seg = "000";

    while ($seg ne $segId) {
        print STDERR "Creating $seg out of $segId\n";

        if (! -e "$path/seg$seg.posDB") {
            runCommand("/bin/sh $path/create.sh $seg") and die "Segment $seg failed.\n";
        }

        $seg++;
        $seg = substr("000$seg", -3);
    }
} elsif ($sge) {

    #  Check if we need to submit pieces of the array, or if we can submit the whole thing.
    #
    my @ap;
    my $wholeThing = 0;

    system("mkdir $path/sgeout") if (! -d "$path/sgeout");

    my $seg = "000";
    while ($seg ne $segId) {
        if (-e "$path/seg$seg.posDB") {
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
        if (runCommand("qsub -cwd -j y -o $path/sgeout/seg\\\$TASK_ID.out -t 1-$segId $sge $path/create.sh")) {
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
                if (runCommand("qsub -cwd -j y -o $path/sgeout/seg\\\$TASK_ID.out -t $st-$ed $sge $path/create.sh")) {
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
} else {
    die "HELP!  I don't know how to run jobs!\n";
}

