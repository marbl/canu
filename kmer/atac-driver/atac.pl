#!/usr/bin/env perl 

$| = 1;
use strict;

use vars qw();

use FindBin;
use lib "$FindBin::Bin/util";

require "run.pl";


my $id1            = undef;
my $seq1           = undef;
my $id2            = undef;
my $seq2           = undef;

my $ATACdir        = undef;
my $GENOMEdir      = "default";   #  Location of genome assemblies
my $MERYLdir       = "default";   #  Location of genome mercount databases

my $mersize        = 20; # the mer size
my $minfill        = 20; # the mimimum fill for a reported match.
my $merlimit       = 1;  # unique mers only
my $maxgap         = 0;  # the maximum substitution gap

my $filtername     = "filter-heavychains.so";
my $filteropts     = "-S 100 -J 100000";

my $numSegments    = 2;
my $numThreads     = 4;

my $merylThreads   = 2;

my $merylOnly      = 0;

my $segmentIDtorun = undef;
my $buildOnly      = undef;

my $execHome       = undef;



if (scalar(@ARGV) < 6) {
    print STDERR "usage: $0 [opts]\n";
    print STDERR "\n";
    print STDERR "    -dir run-directory\n"; # MANDATORY
    print STDERR "\n";
    print STDERR "Sequence specification:  If -seq is supplied, then that\n";
    print STDERR "sequence file is used with the id given by -id.  If there is\n";
    print STDERR "a conflict with an established id, the program exits.\n";
    print STDERR "\n";
    print STDERR "    -id1  id1\n";
    print STDERR "    -seq1 seq1.fasta\n";
    print STDERR "    -id2  id2\n";
    print STDERR "    -seq2 seq2.fasta\n";
    print STDERR "\n";
    print STDERR "NOTE:  id1 is the table, id2 is the stream.\n";
    print STDERR "\n";
    print STDERR "Paths should be FULL PATHS, not relative paths.\n";
    print STDERR "\n";
    print STDERR "    -genomedir path        -- path to the GENOMES directory\n"; # MANDATORY
    print STDERR "    -meryldir  path        -- path to the MERYL directory\n"; # MANDATORY
    print STDERR "    -bindir  path          -- path to the binaries (hack!)\n"; # MANDATORY
    print STDERR "\n";
    print STDERR "    -numsegments s         -- number of segments to do the search in\n";
    print STDERR "    -numthreads t          -- number of threads to use per search\n";
    print STDERR "\n";
    print STDERR "    -merylonly             -- only run the meryl components\n";
    print STDERR "    -merylthreads t        -- number of threads to use for meryl\n";
    print STDERR "\n";
    print STDERR "\n";
    print STDERR "    -samespecies           -- use magic values for same species\n";
    print STDERR "    -crossspecies          -- use guesses for different species\n";
    print STDERR "\n";
    print STDERR "\n";
    print STDERR "ADVANCED OPTIONS:\n";
    print STDERR "\n";
    print STDERR "    -segmentid x           -- only run segment with id x\n";
    print STDERR "                              (don't use unless you really know what it does)\n";
    print STDERR "    -buildonly             -- if defined, don't run the search if the table exists\n";
    print STDERR "\n";
    exit(1);
}


while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if      ($arg eq "-dir") {
        $ATACdir = shift @ARGV;
    } elsif ($arg eq "-id1") {
        $id1 = shift @ARGV;
    } elsif ($arg eq "-seq1") {
        $seq1 = shift @ARGV;
    } elsif ($arg eq "-id2") {
        $id2 = shift @ARGV;
    } elsif ($arg eq "-seq2") {
        $seq2 = shift @ARGV;
    } elsif ($arg eq "-genomedir") {
        $GENOMEdir = shift @ARGV;
    } elsif ($arg eq "-meryldir") {
        $MERYLdir = shift @ARGV;
    } elsif ($arg eq "-numsegments") {
        $numSegments = shift @ARGV;
    } elsif ($arg eq "-numthreads") {
        $numThreads = shift @ARGV;
    } elsif ($arg eq "-bindir") {
        $execHome = shift @ARGV;
    } elsif ($arg eq "-merylonly") {
        $merylOnly = 1;
    } elsif ($arg eq "-merylthreads") {
        $merylThreads = shift @ARGV;
    } elsif ($arg eq "-segmentid") {
        $segmentIDtorun = shift @ARGV;
    } elsif ($arg eq "-buildonly") {
        $buildOnly = 1;
    } elsif ($arg eq "-samespecies") {
        $mersize   = 20; # the mer size
        $merlimit  = 1;  # unique mers only
        $minfill   = 20; # the mimimum fill for a reported match.
        $maxgap    = 0;  # the maximum substitution gap
    } elsif ($arg eq "-samespecies9") {
        $mersize   = 20; # the mer size
        $merlimit  = 9;  # unique mers only
        $minfill   = 20; # the mimimum fill for a reported match.
        $maxgap    = 0;  # the maximum substitution gap
    } elsif ($arg eq "-crossspecies") {
        $mersize   = 18; # the mer size
        $merlimit  = 9;  # unique mers only
        $minfill   = 18; # the mimimum fill for a reported match.
        $maxgap    = 0;  # the maximum substitution gap
    } elsif($arg eq "-filtername") {
        $filtername = shift @ARGV;
    } elsif($arg eq "-filteropts") {
        $filteropts = shift @ARGV;
    } elsif ($arg eq "-mersize") {
        $mersize   = shift @ARGV;
        $minfill   = $mersize;
    } elsif ($arg eq "-merlimit") {
        $merlimit  = shift @ARGV;
    } elsif ($arg eq "-atac") {
        die "-atac unimplemented.  Sorry.\n";
    } else {
        die "unknown option $arg\n";
    }
}


#  Decide on a path to the executables.  This is probably
#  a hack.
#
my $leaff   = "$execHome/leaff";
my $meryl   = "$execHome/meryl";
my $existDB = "$execHome/existDB";
my $seatac  = "$execHome/seatac";

die "Can't run $leaff\n"   if (! -x $leaff);
die "Can't run $meryl\n"   if (! -x $meryl);
die "Can't run $existDB\n" if (! -x $existDB);
die "Can't run $seatac\n"  if (! -x $seatac);

#  Try to find the filter, check here, $execHome and $execHome/../lib
#
if      (-e $filtername) {
    #  nop!
} elsif (-e "$execHome/$filtername") {
    $filtername = "$execHome/$filtername";
} elsif (-e "$execHome/../lib/$filtername") {
    $filtername = "$execHome/../lib/$filtername";
} else {
    print STDERR "Didn't find $filtername.\n";
}

die "Unset GENOMEdir?'\n" if (! defined($GENOMEdir));
die "Unset MERYLdir?'\n"  if (! defined($MERYLdir));
die "Unset ATACdir?'\n"   if (! defined($ATACdir));

die "Can't find the GENOMEdir '$GENOMEdir'\n" if (! -d $GENOMEdir);

system("mkdir $ATACdir")        if (! -d "$ATACdir");
system("mkdir $MERYLdir")       if (! -d "$MERYLdir");


findSources();

my $mercount1 = countMers($id1, $mersize, $merlimit);
my $mercount2 = countMers($id2, $mersize, $merlimit);

my $matches   = "${id1}vs${id2}.k$mersize.u$merlimit.f$minfill.g$maxgap";

#
#  Find the include or exclude mask
#
if (! -e "$ATACdir/$matches.mask.done") {
    my $minFile="min.$mercount1.$mercount2";
    if (! -e "$ATACdir/$minFile.mcdat") {
        print STDERR "Finding the min count between $mercount1 and $mercount2.\n";
        
        my $cmd;
        $cmd  = "$meryl ";
        $cmd .= "-M min ";
        $cmd .= "-s $MERYLdir/$mercount1 ";
        $cmd .= "-s $MERYLdir/$mercount2 ";
        $cmd .= "-o $ATACdir/$minFile ";
        $cmd .= "-stats $ATACdir/$minFile.stats";

        if (runCommand($cmd)) {
            unlink "$ATACdir/$minFile.mcidx";
            unlink "$ATACdir/$minFile.mcdat";
            die "Failed to find the min count between $mercount1 and $mercount2\n";
        }
    }

    die "Failed to make the mask?\n" if (! -e "$ATACdir/$minFile.mcdat");

    #  Decide if we want to use an include mask, or an exclude mask, based
    #  on the estimated size of each.
    #
    #  An include mask is just the 'min' mers found above, while an exclude
    #  mask is 'id1-min' mers.
    #
    my $includeSize = (-s "$ATACdir/$minFile.mcdat");
    my $excludeSize = (-s "$MERYLdir/$mercount1.mcdat") - (-s "$ATACdir/$minFile.mcdat");

    print STDERR "includeSize is about $includeSize\n";
    print STDERR "excludeSize is about $excludeSize\n";

    if ($includeSize < $excludeSize) {
        rename "$ATACdir/$minFile.mcidx", "$ATACdir/$matches.include.mcidx";
        rename "$ATACdir/$minFile.mcdat", "$ATACdir/$matches.include.mcdat";
    } else {
        if (! -e "$ATACdir/$matches.exclude.mcdat") {
            print STDERR "Finding 'exclude' mers!\n";

            my $cmd;
            $cmd  = "$meryl ";
            $cmd .= "-M xor ";
            $cmd .= "-s $MERYLdir/$id2 ";
            $cmd .= "-s $ATACdir/$minFile ";
            $cmd .= "-o $ATACdir/$matches.exclude ";
            $cmd .= "-stats $ATACdir/$matches.exclude.stats";

            if (runCommand($cmd)) {
                unlink "$ATACdir/$matches.exclude.mcidx";
                unlink "$ATACdir/$matches.exclude.mcdat";
                die "Failed to make exclude mers!\n";
            }
        }

        if (-e "$ATACdir/$matches.exclude.mcdat") {
            unlink "$ATACdir/$minFile.mcdat";
            unlink "$ATACdir/$minFile.mcidx";
        } else {
            die "Failed to find exclude mers?\n";
        }
    }

    #  Success!
    #
    open(F, "> $ATACdir/$matches.mask.done");
    close(F);
}

exit(0) if ($merylOnly == 1);

#
#  This is the segmented search routine.  By default, it will segment into two pieces.
#

my $segmentID   = "000";
my @segmentIDs;

open(F, "$leaff -F $MERYLdir/$id1.fasta --partition $numSegments |");
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

#
#  Now, for each segment that hasn't run, run it.
#


#  First seatac run, to build the search tables
#
foreach my $segmentID (@segmentIDs) {
    next if (defined($segmentIDtorun) && ($segmentID ne $segmentIDtorun));

    #  We only want to build in this stage.
    #
    next if (-e "$ATACdir/$matches-segment-$segmentID.table");

    if (! -e "$ATACdir/$matches-segment-$segmentID.build.stats") {
        my $cmd = "";

        $cmd  = "$seatac \\\n";
        $cmd .= "-verbose \\\n";
        $cmd .= "-mersize     $mersize \\\n";
        $cmd .= "-minlength   $minfill \\\n";
        $cmd .= "-maxgap      $maxgap \\\n";
        $cmd .= "-numthreads  $numThreads \\\n";
        $cmd .= "-table       $MERYLdir/$id1.fasta \\\n";
        $cmd .= "-stream      $MERYLdir/$id2.fasta \\\n";
        $cmd .= "-only        $ATACdir/$matches.include \\\n" if (-e "$ATACdir/$matches.include.mcdat");
        $cmd .= "-mask        $ATACdir/$matches.exclude \\\n" if (-e "$ATACdir/$matches.exclude.mcdat");
        $cmd .= "-use         $ATACdir/$matches-segment-$segmentID \\\n";
        $cmd .= "-output      $ATACdir/$matches-segment-$segmentID.matches \\\n";
        $cmd .= "-stats       $ATACdir/$matches-segment-$segmentID.build.stats \\\n";
        $cmd .= "-buildtables $ATACdir/$matches-segment-$segmentID.table \\\n";
        $cmd .= "-filtername  $filtername \\\n" if (defined($filtername));
        $cmd .= "-filteropts  \"-1 $id1 -2 $id2 $filteropts\" ";

        open(F, "> $ATACdir/$matches-$segmentID.build.cmd");
        print F "$cmd\n";
        close(F);

        if (runCommand($cmd)) {
            unlink "$ATACdir/$matches-segment-$segmentID.matches.crash";
            rename "$ATACdir/$matches-segment-$segmentID.matches", "$ATACdir/$matches-segment-$segmentID.matches.crash";
            unlink "$ATACdir/$matches-segment-$segmentID.build.stats.crash";
            rename "$ATACdir/$matches-segment-$segmentID.build.stats", "$ATACdir/$matches-segment-$segmentID.build.stats.crash";
            die "Failed to build tables for $matches-$segmentID\n";
        }
    }
}


#  End early if we are only building
#
if (defined($buildOnly)) {
    print STDERR "Terminating execution because we only wanted to build tables.\n";
    exit(0);
}


#  Check that all the tables are here
#
my $tableNotFound = 0;
foreach my $segmentID (@segmentIDs) {
    if (! -e "$ATACdir/$matches-segment-$segmentID.table") {
        print STDERR "Didn't find a table for $segmentID\n";
        $tableNotFound++;
    }
}
die if ($tableNotFound);


#  Second seatac run, to do the search
#
foreach my $segmentID (@segmentIDs) {
    next if (defined($segmentIDtorun) && ($segmentID ne $segmentIDtorun));

    if (! -e "$ATACdir/$matches-segment-$segmentID.stats") {
        my $cmd = "";

        $cmd  = "$seatac \\\n";
        $cmd .= "-verbose \\\n";
        $cmd .= "-mersize     $mersize \\\n";
        $cmd .= "-minlength   $minfill \\\n";
        $cmd .= "-maxgap      $maxgap \\\n";
        $cmd .= "-numthreads  $numThreads \\\n";
        $cmd .= "-table       $MERYLdir/$id1.fasta \\\n";
        $cmd .= "-stream      $MERYLdir/$id2.fasta \\\n";
        $cmd .= "-only        $ATACdir/$matches.include \\\n" if (-e "$ATACdir/$matches.include.mcdat");
        $cmd .= "-mask        $ATACdir/$matches.exclude \\\n" if (-e "$ATACdir/$matches.exclude.mcdat");
        $cmd .= "-use         $ATACdir/$matches-segment-$segmentID \\\n";
        $cmd .= "-output      $ATACdir/$matches-segment-$segmentID.matches \\\n";
        $cmd .= "-usetables   $ATACdir/$matches-segment-$segmentID.table \\\n";
        $cmd .= "-stats       $ATACdir/$matches-segment-$segmentID.stats \\\n";
        $cmd .= "-filtername  $filtername \\\n" if (defined($filtername));
        $cmd .= "-filteropts  \"-1 $id1 -2 $id2 $filteropts\" ";

        #  Prevent me from overwriting a run in progress
        #
        if (-e "$ATACdir/$matches-segment-$segmentID.matches") {
            die "WARNING:  Matches already exist!  Exiting!\n";
        }

        open(F, "> $ATACdir/$matches-$segmentID.cmd");
        print F "$cmd\n";
        close(F);

        if (runCommand($cmd)) {
            unlink "$ATACdir/$matches-segment-$segmentID.matches.crash";
            rename "$ATACdir/$matches-segment-$segmentID.matches", "$ATACdir/$matches-segment-$segmentID.matches.crash";
            unlink "$ATACdir/$matches-segment-$segmentID.stats.crash";
            rename "$ATACdir/$matches-segment-$segmentID.stats", "$ATACdir/$matches-segment-$segmentID.stats.crash";
            die "Failed to run $matches-$segmentID\n";
        }
    }
}

#  End early if the segment id to run is defined.
#
if (defined($segmentIDtorun)) {
    print STDERR "Terminating execution because a specific segmentID was supplied.\n";
    exit(0);
}



#
#  Join and sort the matches
#
if (! -e "$ATACdir/$matches.matches.sorted") {
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

    open(ATACFILE, "> $ATACdir/$matches.matches.sorted");
    print ATACFILE "!format atac 1.0\n";
    print ATACFILE  "#\n";
    print ATACFILE  "# Legend:\n";
    print ATACFILE  "#\n";
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
    print ATACFILE  "#\n";
    print ATACFILE "/assemblyId1=$id1\n";
    print ATACFILE "/assemblyId2=$id2\n";
    print ATACFILE "/assemblyFile1=$MERYLdir/$id1.fasta\n";
    print ATACFILE "/assemblyFile2=$MERYLdir/$id2.fasta\n";

    #  We used to trim off the fasta from the filename...why?
    my $seq1trimmed = $seq1;
    my $seq2trimmed = $seq2;
    $seq1trimmed = $1 if ($seq1trimmed =~ m/(.*).fasta$/);
    $seq2trimmed = $1 if ($seq2trimmed =~ m/(.*).fasta$/);

    print ATACFILE "/rawMatchMerSize=$mersize\n";
    print ATACFILE "/rawMatchMerMaxDegeneracy=$merlimit\n";
    print ATACFILE "/rawMatchAllowedSubstutionBlockSize=$maxgap\n";
    print ATACFILE "/rawMatchMinFillSize=$minfill\n";

    print ATACFILE "/heavyChainsOn=1\n";
    print ATACFILE "/heavyMaxJump=100000\n";
    print ATACFILE "/heavyMinFill=100\n";

    print ATACFILE "/matchExtenderOn=1\n";

    print ATACFILE "/uniqueFilterOn=1\n";
    print ATACFILE "/fillIntraRunGapsOn=1\n";

    if (0){
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

    #print ATACFILE "/matchesFile=$ATACdir/$matches.matches.sorted.bz2\n";

    close(ATACFILE);




    #  The original sort used keys seq1, seq2, pos1, pos2.  The next
    #  consumer of this data, matchextender certainly benefits from
    #  the first two keys, and probably needs the last two....nope, it
    #  does an internal sort given all matches for a pair of
    #  sequences.

    system("mkdir ${ATACdir}/tmp")  if (! -d "${ATACdir}/tmp");

    if (runCommand("cat $mfiles | grep ^M | sort -y -T $ATACdir/tmp -k 5,5 -k 9,9 >> $ATACdir/$matches.matches.sorted")) {
        die "Failed to sort $ATACdir!\n";
    }

    system("rm -rf ${ATACdir}/tmp");
}


#  Run matchExtender
#
if ((! -e "$ATACdir/$matches.matches.sorted.extended") ||
    (0 == -s "$ATACdir/$matches.matches.sorted.extended")) {
    my $cmd;
    $cmd  = "$execHome/MatchExtender -v ";
    $cmd .= "$ATACdir/$matches.matches.sorted ";
    $cmd .= "$ATACdir/$matches.matches.sorted.extended ";

    if (runCommand($cmd)) {
        print STDERR "MatchExtender failed.\n";
        exit(1);
    }
}










#  Read the nickname file, set up symlinks to the data sources
#
sub findSources {
    my %GENOMEaliases;

    #  Read all the *.atai files in the genome directory, save only
    #  those nicknames that have actual files associated with them.
    #  This lets us have multiple collections of assemblies, and also
    #  lets us move the directory around (e.g., for running on a
    #  laptop).
    #
    open(F, "cat $GENOMEdir/*.atai |") or die "Can't cat $GENOMEdir/*.atai\n";
    while (<F>) {
        chomp;

        if (m/^!\s*format\s+atac\s+(.*)$/) {
            print STDERR "Found format $1\n";
        } elsif (m/^S\s+(\S+)\s+(\S+)$/) {
            $GENOMEaliases{$1} = $2 if (-e $2);
        } else {
            #die "Error parsing genome description.\n  '$_'\n";
        }
    }
    close(F);

    #  If the user gave both an id and a sequence, make sure that
    #  the id is distinct.
    #
    die "No id1 supplied!\n" if (!defined($id1));
    die "No id2 supplied!\n" if (!defined($id2));

    die "id1 = '$id1' is already used by sequence '$GENOMEaliases{$id1}'\n" if (defined($GENOMEaliases{$id1}) && defined($seq1));
    die "id2 = '$id2' is already used by sequence '$GENOMEaliases{$id2}'\n" if (defined($GENOMEaliases{$id2}) && defined($seq2));

    $seq1 = $GENOMEaliases{$id1} if (!defined($seq1));
    $seq2 = $GENOMEaliases{$id2} if (!defined($seq2));

    die "Unknown alias $id1.\n" if (!defined($seq1));
    die "Unknown alias $id2.\n" if (!defined($seq2));

    die "File '$seq1' doesn't exist for alias $id1.\n" if (! -e $seq1);
    die "File '$seq2' doesn't exist for alias $id2.\n" if (! -e $seq2);
    
    system("ln -s $seq1 $MERYLdir/$id1.fasta") if (! -e "$MERYLdir/$id1.fasta");
    system("ln -s $seq2 $MERYLdir/$id2.fasta") if (! -e "$MERYLdir/$id2.fasta");

    system("ln -s ${seq1}idx $MERYLdir/$id1.fastaidx") if (! -e "$MERYLdir/$id1.fastaidx") && (-e "${seq1}idx");
    system("ln -s ${seq2}idx $MERYLdir/$id2.fastaidx") if (! -e "$MERYLdir/$id2.fastaidx") && (-e "${seq2}idx");
}


#  Check that meryl is finished for each of the inputs
#
sub countMers {
    my ($id, $mersize, $merlimit) = @_;

    #  Using "-H 32" is needed if the two sequences aren't about the
    #  same order of magnitude in size.  This value is appropriate for
    #  sequences that are genome size.

    if (! -e "$MERYLdir/$id.mcdat") {
        my $cmd;
        $cmd  = "$meryl -B -C ";
        $cmd .= "-threads $merylThreads ";
        $cmd .= "-m $mersize ";
        $cmd .= "-s $MERYLdir/$id.fasta ";
        $cmd .= "-o $MERYLdir/$id ";
        $cmd .= "-stats $MERYLdir/$id.stats";
        if (runCommand($cmd)) {
            unlink "$MERYLdir/$id.mcidx";
            unlink "$MERYLdir/$id.mcdat";
            die "Failed to count mers in $id\n";
        }
    }

    if (! -e "$MERYLdir/$id.le$merlimit.mcdat") {
        my $cmd;
        $cmd  = "$meryl -v ";
        $cmd .= "-M lessthanorequal $merlimit ";
        $cmd .= "-s $MERYLdir/$id ";
        $cmd .= "-o $MERYLdir/$id.le$merlimit ";
        $cmd .= "-stats $MERYLdir/$id.le$merlimit.stats";
        if (runCommand($cmd)) {
            unlink "$MERYLdir/$id.le$merlimit.mcidx";
            unlink "$MERYLdir/$id.le$merlimit.mcdat";
            die "Failed to count mers lessthanorequal $merlimit in $id\n";
        }
    }

    return "$id.le$merlimit";
}

