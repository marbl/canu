#!/usr/bin/perl

#  Runs snapper2 on SGE, splitting both the genome and query sequences.

use FindBin;

my $genome = "";
my $query  = "";
my $dir    = "";
my $mask   = 1000;
my $gseg   = 32;
my $qseg   = 32;
my $check  = undef;

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
    } elsif ($arg =~ m/^-gseg/) {
        $gseg   = shift @ARGV;
    } elsif ($arg =~ m/^-qseg/) {
        $qseg   = shift @ARGV;
    } elsif ($arg =~ m/^-check/) {
        $check = 1;
    } else {
        print STDERR "unknown option '$arg'\n";
    }
}

#  If we're checking the results, assume we are in the correct dir,
#  that the gen and qry dirs exist.
#
if (defined($check)) {
    my ($gen, $qry) = countSequences($dir);

    for (my $g=1; $g<=$gen; $g++) {
        for (my $q=1; $q<=$qry; $q++) {
            $g = substr("000$g", -3);
            $q = substr("000$q", -3);

            if (0) {
                print STDERR "Why am I doing this?  I'm supposed to be checking overlap, not snapper....\n";
                exit(1);
            }
        }
    }
}




if (!defined($genome) || !defined($query) || !defined($dir)) {
    print STDERR "usage: $0 [arg]\n";
    print STDERR "  -genome x.fasta\n";
    print STDERR "  -query  q.fasta\n";
    print STDERR "  -dir    /path/to/work\n";
    print STDERR "  -mask   kmer-limit           (def: 1000)\n";
    print STDERR "  -gseg   gseg                 (def: 16 segs, see leaff for format\n";
    print STDERR "  -qseg   qseg                 (def: 16 segs, see leaff for format\n";
    print STDERR "  -check                       (check a run, assume we are in the /path/to/work\n";
    exit(1);
}

die "Can't find genome '$genome'\n" if (! -e $genome);
die "Can't find queries '$query'\n"  if (! -e $query);

system("mkdir $dir")      if (! -d $dir);
die "Can't find '$dir'\n" if (! -d $dir);

if (! -e "$dir/gen/gen.partitioned") {
    system("mkdir $dir/gen")  if (! -d "$dir/gen");
    system("$bin/leaff -F $genome --partition $dir/gen/gen $gseg");
    open(F, "> $dir/gen/gen.partitioned");
    close(F);
}

if (! -e "$dir/qry/qry.partitioned") {
    system("mkdir $dir/qry")  if (! -d "$dir/qry");
    system("$bin/leaff -F $query  --partition $dir/qry/qry $qseg");
    open(F, "> $dir/qry/qry.partitioned");
    close(F);
}

#  Build indexes for everyone -- this prevents the grid jobs from
#  racing to build them.  And it lets us count how many jobs to
#  submit.
#
my ($gen, $qry) = countSequences($dir);

open(F, "> $dir/run.sh");
print F "#!/bin/sh\n";
print F "PIECE=`expr \$SGE_TASK_ID - 1`\n";
print F "GPIECE=`expr \$PIECE % $gen + 1`\n";
print F "QPIECE=`expr \$PIECE / $gen + 1`\n";
print F "GPIECE=`printf %03d \$GPIECE`\n";
print F "QPIECE=`printf %03d \$QPIECE`\n";
print F "scratchname=/scratch/\$\$-\$GPIECE-\$QPIECE\n";
print F "\n";
print F "ulimit -c 0\n";
print F "#rm /scratch/[0-9]*-[0-9]*-[0-9]*\n";
print F "#echo $GPIECE $QPIECE $PIECE\n";
print F "\n";
print F "if [ -e $dir/map-gen$GPIECE-qlt$QPIECE.success ] ; then\n";
print F "  echo map-gen$GPIECE-qlt$QPIECE already done\n";
print F "  exit\n";
print F "fi\n";
print F "\n";
print F "$bin/snapper2 -verbose \\\n";
print F "  -mersize 22 -merskip 0 \\\n";
print F "  -minhitlength 22 -minhitcoverage 0.0 \\\n";
#print F "  -setfilter 0.1500 0.1500 0.2500 \\\n";
print F "  -validate $dir/map-gen\$GPIECE-qlt\$QPIECE.validate \\\n";
print F "  -genomic $dir/gen/gen-\$GPIECE.fasta \\\n";
print F "  -queries $dir/qry/qry-\$QPIECE.fasta \\\n";
print F "  -ignore  $mask \\\n";
#print F "  -output  \$scratchname \\\n";
print F "  -noaligns \\\n";
print F "  -numthreads 2 \\\n";
print F "  -minmatchidentity 90 \\\n";
print F "  -minmatchcoverage 4 \\\n";
print F "  -loaderhighwatermark 1024 \\\n";
print F "| \\\n";
print F "bzip2 -9vc > \$scratchname.bz2 \\\n";
print F "&& \\\n";
print F "mv \$scratchname.bz2 $dir/map-gen\$GPIECE-qlt\$QPIECE.sim4db.bz2 \\\n";
print F "&& \\\n";
print F "touch $dir/map-gen\$GPIECE-qlt\$QPIECE.success\n";
close(F);

my $numJobs = $gen * $qry;

print STDOUT "qsub -pe thread 2 -t 1-$numJobs -p -50 -j y -o $dir/map-\\\$TASK_ID $dir/run.sh\n";





sub countSequences {
    my $dir = shift @_;
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

    return($gen, $qry);
}
