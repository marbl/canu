#!/usr/local/bin/perl -w

use strict;
use Config;  #  for @signame

#  Usage:  atac.pl run-directory id1 id2
#
if (scalar(@ARGV) != 3) {
    print STDERR "usage: $0 run-directory id1 id2\n";
    exit(1);
}

my $ATACdir = shift @ARGV;
my $id1     = shift @ARGV;
my $id2     = shift @ARGV;

my $GENOMEdir     = "/prod/IR05/GENOMES";          #  Location of genome assemblies
my $MERYLdir      = "/prod/IR08/walenz/hg5/data";  #  Location of genome mercount databases
my %GENOMEaliases;

my $leaff   = "/work/assembly/walenzbp/releases/leaff";
my $meryl   = "/work/assembly/walenzbp/releases/meryl";
my $existDB = "/work/assembly/walenzbp/releases/existDB";
my $sge     = "/work/assembly/walenzbp/releases/searchGENOMEexactly";

die "Can't run $leaff\n"   if (! -x $leaff);
die "Can't run $meryl\n"   if (! -x $meryl);
die "Can't run $existDB\n" if (! -x $existDB);
die "Can't run $sge\n"     if (! -x $sge);

die "Can't find the GENOMEdir '$GENOMEdir'\n" if (! -d $GENOMEdir);
die "Can't find the assembly descriptions '$GENOMEdir/assemblies.atai'\n" if (! -e "$GENOMEdir/assemblies.atai");
die "Can't find the MERYLdir '$MERYLdir'\n" if (! -d $MERYLdir);

system("mkdir $ATACdir") if (! -d $ATACdir);


my $mersize   = 20; # the mer size
my $merlimit  = 1;  # unique mers only
my $minfill   = 20; # the mimimum fill for a reported match.
my $maxgap    = 0;  # the maximum substitution gap

my $numSegments = 2;

findSources($id1, $id2);

my $mercount1 = countMers($id1, $mersize, $merlimit);
my $mercount2 = countMers($id2, $mersize, $merlimit);

my $matches   = "$id1-vs-$id2.k$mersize.u$merlimit.f$minfill.g$maxgap";

#
#  Find the include or exclude mask
#
if (! -e "$ATACdir/.mask.done") {

    if (! -e "$ATACdir/min.$mercount1.$mercount2.mcdat") {
        print STDERR "Finding the min count between $mercount1 and $mercount2.\n";
        if (runCommand("$meryl -v -M min -s $MERYLdir/$mercount1 -s $MERYLdir/$mercount2 -o $ATACdir/min.$mercount1.$mercount2")) {
            unlink "$ATACdir/min.$mercount1.$mercount2.mcidx";
            unlink "$ATACdir/min.$mercount1.$mercount2.mcdat";
            die "Failed to find the min count between $mercount1 and $mercount2\n";
        }
    }

    die "Failed to make the mask?\n" if (! -e "$ATACdir/min.$mercount1.$mercount2.mcdat");

    #  Decide if we want to use an include mask, or an exclude mask, based
    #  on the estimated size of each.
    #
    #  An include mask is just the 'min' mers found above, while an exclude
    #  mask is 'id1-min' mers.
    #
    my $includeSize = (-s "$ATACdir/min.$mercount1.$mercount2.mcdat");
    my $excludeSize = (-s "$MERYLdir/$mercount1.mcdat") - (-s "$ATACdir/min.$mercount1.$mercount2.mcdat");

    print STDERR "includeSize is about $includeSize\n";
    print STDERR "excludeSize is about $excludeSize\n";

    #
    #  Since we usually run multiple copies of the search, and since building
    #  the existDB structure takes > thirty minutes, we pre-build it.
    #

    if ($includeSize < $excludeSize) {

        if (! -e "$ATACdir/$matches.include.existDB") {
            print STDERR "Building 'include' existDB structure.\n";
            if (runCommand("$existDB -m 20 -t 19 $ATACdir/min.$mercount1.$mercount2 $ATACdir/$matches.include.existDB")) {
                unlink "$ATACdir/$matches.include.existDB";
                die "Failed to make include existDB?\n";
            }
        }

        die "Failed to make include existDB?\n" if (! -e "$ATACdir/$matches.include.existDB");
    } else {

        if (! -e "$ATACdir/$matches.exclude.mcdat") {
            print STDERR "Finding 'exclude' mers!\n";
            if (runCommand("./meryl -v -M xor -s $MERYLdir/$id1 -s $ATACdir/min.$mercount1.$mercount2 -o $ATACdir/$matches.exclude")) {
                unlink "$ATACdir/$matches.exclude.mcidx";
                unlink "$ATACdir/$matches.exclude.mcdat";
                die "Failed to make exclude mers!\n";
            }
        }

        die "Failed to find exclude mers?\n" if (! -e "$ATACdir/$matches.exclude.mcdat");

        if (! -e "$ATACdir/$matches.exclude.existDB") {
            print STDERR "Building 'exclude' existDB structure.\n";
            if (runCommand("$existDB -m 20 -t 19 $ATACdir/$matches.exclude $ATACdir/$matches.exclude.existDB")) {
                unlink "$ATACdir/$matches.exclude.existDB";
                die "Failed to make exclude existDB?\n";
            }
        }

        die "Failed to make exclude existDB?\n" if (! -e "$ATACdir/$matches.exclude.existDB");
    }

    #  Success!
    #
    open(F, "> $ATACdir/.mask.done");
    close(F);

    #  Clean up the temporary files
    #
    #unlink "$ATACdir/min.$mercount1.$mercount2.mcdat";
    #unlink "$ATACdir/min.$mercount1.$mercount2.mcidx";
    #unlink "$ATACdir/$matches.exclude.mcdat";
    #unlink "$ATACdir/$matches.exclude.mcidx";
}


exit(0);


#  This is the segmented search routine.  By default, it will segment into two pieces.
#
#  $id1 is used as the "query" sequences
#  $id2 is used for the table
#

my $segmentID   = "000";
my @segmentIDs;

open(F, "$leaff -F $MERYLdir/$id2.fasta --partition $numSegments |");
$numSegments = <F>;
while(<F>) {
    my $segments = "";
    my @pieces   = split '\s+', $_;
    my $memory   = shift @pieces;

    foreach my $piece (@pieces) {
        if ($piece =~ m/(\d+)\(\d+\)/) {
            $segments .= "$1\n";
        } else {
            die "Error parsing segment: $piece\n";
        }
    }

    open(S, "> $ATACdir/$matches-segment-$segmentID");
    print S $segments;
    close(S);

    push @segmentIDs, $segmentID;

    $segmentID++;
}
close(F);



#  Now, for each segment that hasn't run, run it.
#
foreach my $segmentID (@segmentIDs) {
    if (! -e "$ATACdir/$matches-$segmentID.stats") {
        my $cmd = "";

        $cmd  = "$sge -verbose ";
        $cmd .= "-loaderhighwatermark 2 ";
        $cmd .= "-mersize $mersize ";
        $cmd .= "-singlelength $minfill ";
        $cmd .= "-maxgap $maxgap ";
        $cmd .= "-numthreads 4 ";
        $cmd .= "-genomic $MERYLdir/$id2.fasta ";
        $cmd .= "-cdna $MERYLdir/$id1.fasta ";
        $cmd .= "-only $ATACdir/$matches.include.existDB " if (-e "$ATACdir/$matches.include.existDB");
        $cmd .= "-mask $ATACdir/$matches.exclude.existDB " if (-e "$ATACdir/$matches.exclude.existDB");
        $cmd .= "-use $ATACdir/$matches-segment-$segmentID ";
        $cmd .= "-output $ATACdir/$matches-segment-$segmentID.matches ";
        $cmd .= "-stats $ATACdir/$matches-segment-$segmentID.stats ";

        open(F, "> $ATACdir/$matches-$segmentID.cmd");
        print F "$cmd\n";
        close(F);

        if (runCommand($cmd)) {
            unlink "$ATACdir/$matches-segment-$segmentID.matches";
            unlink "$ATACdir/$matches-segment-$segmentID.stats";
            die "Failed to run $matches-$segmentID\n";
        }
    }
}

#
#  Join and sort the matches
#
if (! -e "$ATACdir/$matches.matches") {
    my $mfiles;

    #  Check that each search finished, and build a list of all the match files.
    #
    foreach my $segmentID (@segmentIDs) {
        if (-e "$ATACdir/$matches-segment-$segmentID.stats") {
            $mfiles .= "$ATACdir/$matches-segment-$segmentID.matches ";
        } else {
            die "$ATACdir/$matches-segment-$segmentID.matches failed to complete.\n";
        }
    }

    if (runCommand("cat $mfiles | sort -y -T $ATACdir/sortingjunk -k 3n -k 7n > $ATACdir/$matches.matches.sorted &")) {
        die "Failed to sort $ATACdir!\n";
    }

    #system("rm -f $mfiles");
}


if (! -e "$ATACdir/${id1}vs${id2}-C13.atac") {
    open(ATACFILE, "> $ATACdir/${id1}vs${id2}-C13.atac");
    print ATACFILE  "!format atac 1.0\n";
    print ATACFILE  "# Legend:\n";
    print ATACFILE  "# Field 0: the row class\n";
    print ATACFILE  "# Field 1: the match type u=ungapped, x=exact, ....\n";
    print ATACFILE  "# Field 2: the match instance index\n";
    print ATACFILE  "# Field 3: the parent index\n";
    print ATACFILE  "# Field 4: the FASTA sequence id in the first assembly\n";
    print ATACFILE  "# Field 5: the offset from the start of the sequence for the match\n";
    print ATACFILE  "# Field 6: the length of the match in the first assembly\n";
    print ATACFILE  "# Field 7: the orientation of the match sequence in the first assembly.\n";
    print ATACFILE  "# Field 8: the FASTA sequence id for the second assembly\n";
    print ATACFILE  "# Field 9: the offset from the start of the sequence for the match\n";
    print ATACFILE  "# Field 10: the length of the match in the second assembly\n";
    print ATACFILE  "# Field 11: the orientation of the match sequence in the second assembly.\n";
    print ATACFILE "/assemblyFilePrefix1=$GENOMEaliases{$id1}\n";
    print ATACFILE "/assemblyFilePrefix2=$GENOMEaliases{$id2}\n";
    print ATACFILE "/assemblyId1=$id1\n";
    print ATACFILE "/assemblyId2=$id2\n";
    print ATACFILE "/rawMatchMerSize=$mersize\n";
    print ATACFILE "/rawMatchMerMaxDegeneracy=$merlimit\n";
    print ATACFILE "/rawMatchAllowedSubstutionBlockSize=$maxgap\n";
    print ATACFILE "/rawMatchMinFillSize=$minfill\n";

    #/chain-greedy-on=1
    #/chain-consv-on=1
    #/chain-global-on=1

    print ATACFILE "/heavyChainsOn=1\n";

    #print ATACFILE "/heavyMaxJump=100000\n";
    #print ATACFILE "/heavyMinFill=100\n";

    print ATACFILE "/matchExtenderOn=1\n";

    if(0){
        # The non-default parameters for Mouse versus Rat.
        print ATACFILE "/matchExtenderMaxMMBlock=5\n";
        print ATACFILE "/matchExtenderMaxNbrPathMM=25\n";
        print ATACFILE "/matchExtenderMaxNbrSep=100\n";
        print ATACFILE "/matchExtenderMinBlockSep=5\n";
        print ATACFILE "/matchExtenderMinEndRunLen=4\n";
        print ATACFILE "/matchExtenderMinIdentity=0.7\n";
        print ATACFILE "/globalMatchMinSize=20\n";
        print ATACFILE "/fillIntraRunGapsErate=0.30\n";
    }

    print ATACFILE "/uniqueFilterOn=1\n";
    print ATACFILE "/fillIntraRunGapsOn=1\n";
    print ATACFILE "/matchesFile=$ATACdir/$matches.matches.sorted\n";
    close(ATACFILE);
}




#  Utility to run a command and check the exit status
#
sub runCommand {
    my $cmd = shift @_;
    my $rc = 0xffff & system($cmd);

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!
    #
    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    }

    my $error = "ERROR: $cmd\n        failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } elsif ($rc > 0x80) {
        $rc >>= 8;
        $error .= "exit status $rc\n";
    } else {
        if ($rc & 0x80) {
            $rc &= ~0x80;
            $error .= "coredump from ";
        }
        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc]\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}


#  Read the nickname file, set up symlinks to the data sources
#
sub findSources {
    my $id1 = shift @_;
    my $id2 = shift @_;

    #  Read the assemblies.atai file to generate a mapping of datasource and nickname.
    #
    open(F, "< $GENOMEdir/assemblies.atai") or die "Can't file $GENOMEdir/assemblies.atai\n";
    while (<F>) {
        chomp;

        if (m/^S\s+(\S+)\s+(\S+)$/) {
            $GENOMEaliases{$1} = $2;
        } else {
            die "Error parsing assemblies.atai.\n  '$_'\n";
        }
    }
    close(F);

    my $g1 = $GENOMEaliases{$id1};
    my $g2 = $GENOMEaliases{$id2};

    die "Unknown alias $id1.\n" if (!defined($g1));
    die "Unknown alias $id2.\n" if (!defined($g2));

    die "File '$g1' doesn't exist for alias $id1.\n" if (! -e $g1);
    die "File '$g2' doesn't exist for alias $id2.\n" if (! -e $g2);
    
    system("ln -s $g1 $MERYLdir/$id1.fasta") if (! -e "$MERYLdir/$id1.fasta");
    system("ln -s $g2 $MERYLdir/$id2.fasta") if (! -e "$MERYLdir/$id2.fasta");

    system("ln -s ${g1}idx $MERYLdir/$id1.fastaidx") if (! -e "$MERYLdir/$id1.fastaidx") && (-e "${g1}idx");
    system("ln -s ${g2}idx $MERYLdir/$id2.fastaidx") if (! -e "$MERYLdir/$id2.fastaidx") && (-e "${g2}idx");
}


#  Check that meryl is finished for each of the inputs
#
sub countMers {
    my ($id, $mersize, $merlimit) = @_;

    if (! -e "$MERYLdir/$id.mcdat") {
        if (runCommand("$meryl -v -B -C -m $mersize -t 27 -s $MERYLdir/$id.fasta -o $MERYLdir/$id")) {
            unlink "$MERYLdir/$id.mcidx";
            unlink "$MERYLdir/$id.mcdat";
            die "Failed to count mers in $id\n";
        }
    }

    if (! -e "$MERYLdir/$id.le$merlimit.mcdat") {
        if (runCommand("$meryl -v -M lessthanorequal $merlimit -s $MERYLdir/$id -o $MERYLdir/$id.le$merlimit")) {
            unlink "$MERYLdir/$id.le$merlimit.mcidx";
            unlink "$MERYLdir/$id.le$merlimit.mcdat";
            die "Failed to count mers lessthanorequal $merlimit in $id\n";
        }
    }

    return "$id.le$merlimit";
}

