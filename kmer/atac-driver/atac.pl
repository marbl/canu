#!/usr/bin/env perl
#
# This file is part of A2Amapper.
# Copyright (c) 2005 J. Craig Venter Institute
# Author: Brian Walenz
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

use strict;
use FindBin;

my $id1               = undef;
my $seq1              = undef;
my $id2               = undef;
my $seq2              = undef;

my $ATACdir           = undef;
my $GENOMEdir         = "default";   #  Location of genome assemblies
my $MERYLdir          = "default";   #  Location of genome mercount databases

my $BINdir            = "$FindBin::Bin";
my $LIBdir            = "$FindBin::Bin/../lib";

my $mersize           = 20; # the mer size
my $minfill           = 20; # the mimimum fill for a reported match.
my $merlimit          = 1;  # unique mers only
my $maxgap            = 0;  # the maximum substitution gap

# annotates the resulting atac file with parameters
# for cross species, also sets match extender options
my $crossSpecies      = 0;  

my $matchExtenderOpts = "";

my $filtername        = "$LIBdir/filter-heavychains.so";
my $filteropts        = "-S 100 -J 100000";

my $numSegments       = 2;
my $numThreads        = 4;

my $merylThreads      = 2;

my $merylOnly         = 0;

#  Check that we have everything we need to run
#
my $leaff             = "$BINdir/leaff";
my $meryl             = "$BINdir/meryl";
my $existDB           = "$BINdir/existDB";
my $seatac            = "$BINdir/seatac";
my $chainer           = "$BINdir/AtacDriver.py";

die "Can't run $leaff\n"        if (! -x $leaff);
die "Can't run $meryl\n"        if (! -x $meryl);
die "Can't run $existDB\n"      if (! -x $existDB);
die "Can't run $seatac\n"       if (! -x $seatac);
die "Can't find $chainer\n"     if (! -e $chainer);
die "Can't find $filtername\n"  if (! -e $filtername);



parseArgs();
findSources();

my $mercount1 = countMers($id1, $mersize, $merlimit);
my $mercount2 = countMers($id2, $mersize, $merlimit);

my $matches   = "${id1}vs${id2}.k$mersize.u$merlimit.f$minfill.g$maxgap";

buildMask();

my @segmentIDs = findHits();

sortMatches(@segmentIDs);
extendMatches();

makeChains();
makeClumps();

print STDERR "\n";
print STDERR "Finished!  Output is:\n";
print STDERR "  $ATACdir/$matches.atac\n";
print STDERR "  $ATACdir/$matches.atac.clumps5000\n";


#  Subroutines below!


sub usage {
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
    print STDERR "    -meryldir  path        -- path to the MERYL directory\n";   # MANDATORY
    print STDERR "    -bindir  path          -- path to the binaries (hack!)\n";  # MANDATORY
    print STDERR "    -srcdir  path          -- path to the source   (hack!)\n";  # MANDATORY
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
    print STDERR "\n";
    exit(1);
}

sub parseArgs {

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
        } elsif ($arg eq "-merylonly") {
            $merylOnly = 1;
        } elsif ($arg eq "-merylthreads") {
            $merylThreads = shift @ARGV;
        } elsif ($arg eq "-samespecies") {
            $mersize      = 20; # the mer size
            $merlimit     = 1;  # unique mers only
            $minfill      = 20; # the mimimum fill for a reported match.
            $maxgap       = 0;  # the maximum substitution gap
        } elsif ($arg eq "-samespecies9") {
            $mersize      = 20; # the mer size
            $merlimit     = 9;  # mostly unique mers only
            $minfill      = 20; # the mimimum fill for a reported match.
            $maxgap       = 0;  # the maximum substitution gap
        } elsif ($arg eq "-crossspecies1") {
            $mersize      = 20; # the mer size
            $merlimit     = 1;  # mostly unique mers only
            $minfill      = 20; # the mimimum fill for a reported match.
            $maxgap       = 0;  # the maximum substitution gap
            $crossSpecies = 1;  # extra parameters in the atac file
        } elsif ($arg eq "-crossspecies") {
            $mersize      = 18; # the mer size
            $merlimit     = 9;  # mostly unique mers only
            $minfill      = 18; # the mimimum fill for a reported match.
            $maxgap       = 0;  # the maximum substitution gap
            $crossSpecies = 1;  # extra parameters in the atac file
        } elsif($arg eq "-filtername") {
            $filtername = shift @ARGV;
        } elsif($arg eq "-filteropts") {
            $filteropts = shift @ARGV;
        } elsif ($arg eq "-mersize") {
            $mersize   = shift @ARGV;
            $minfill   = $mersize;
        } elsif ($arg eq "-merlimit") {
            $merlimit  = shift @ARGV;
        } elsif ($arg eq "-justtestingifitworks") {
            exit(0);
        } else {
            die "unknown option $arg\n";
        }
    }

    if (!defined($id1) ||
        !defined($id2)) {
        usage();
    }

    my $pwd = `pwd`;
    $pwd =~ s/^\s+//;
    $pwd =~ s/\s+$//;

    $GENOMEdir = "$pwd/$GENOMEdir" if ($GENOMEdir !~ m!^/!);
    $MERYLdir  = "$pwd/$MERYLdir"  if ($MERYLdir !~ m!^/!);
    $ATACdir   = "$pwd/$ATACdir"   if ($ATACdir !~ m!^/!);

    die "Unset GENOMEdir?'\n" if (! defined($GENOMEdir));
    die "Unset MERYLdir?'\n"  if (! defined($MERYLdir));
    die "Unset ATACdir?'\n"   if (! defined($ATACdir));

    if (!defined($seq1) || (!defined($seq2))) {
        die "Can't find the GENOMEdir '$GENOMEdir'\n" if (! -d $GENOMEdir);
    }
    if (defined($seq1)) {
        $seq1 = "$pwd/$seq1" if ($seq1 !~ m!^/!);
    }
    if (defined($seq2)) {
        $seq2 = "$pwd/$seq2" if ($seq1 !~ m!^/!);
    }

    system("mkdir $ATACdir")        if (! -d "$ATACdir");
    system("mkdir $MERYLdir")       if (! -d "$MERYLdir");
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
    if (-d $GENOMEdir) {
        #  What?  No GENOMEdir?  The main already checked that we know both
        #  sequence files.  Plus, we'd just fail below.

        open(A, "ls $GENOMEdir |");
        while (<A>) {
            chomp;
            if (m/\.atai$/) {
                my $ataifile = "$GENOMEdir/$_";
                open(F, "< $ataifile") or die "Can't open '$ataifile'\n";
                while (<F>) {
                    chomp;

                    if (m/^!\s*format\s+atac\s+(.*)$/) {
                        print STDERR "Found format $1\n";
                    } elsif (m/^S\s+(\S+)\s+(\S+)$/) {
                        if (-e $2) {
                            $GENOMEaliases{$1} = $2;
                        } else {
                            print STDERR "WARNING:  File '$2' not found for alias '$1'.\n";
                        }
                    } else {
                        #die "Error parsing genome description.\n  '$_'\n";
                    }
                }
                close(F);
            }
        }
    }
    close(A);

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

    if (! -e "$MERYLdir/$id.ms$mersize.mcdat") {
        my $cmd;
        $cmd  = "$meryl -B -C ";
        $cmd .= "-threads $merylThreads ";
        $cmd .= "-m $mersize ";
        $cmd .= "-s $MERYLdir/$id.fasta ";
        $cmd .= "-o $MERYLdir/$id.ms$mersize ";
        $cmd .= "-stats $MERYLdir/$id.ms$mersize.stats";
        #die "why rebuild $MERYLdir/$id.ms$mersize.mcdat\n";
        if (runCommand($cmd)) {
            unlink "$MERYLdir/$id.ms$mersize.mcidx";
            unlink "$MERYLdir/$id.ms$mersize.mcdat";
            die "Failed to count mers in $id\n";
        }
    }

    if (! -e "$MERYLdir/$id.ms$mersize.le$merlimit.mcdat") {
        my $cmd;
        $cmd  = "$meryl -v ";
        $cmd .= "-M lessthanorequal $merlimit ";
        $cmd .= "-s $MERYLdir/$id.ms$mersize ";
        $cmd .= "-o $MERYLdir/$id.ms$mersize.le$merlimit ";
        $cmd .= "-stats $MERYLdir/$id.ms$mersize.le$merlimit.stats";
        #die "why rebuild $MERYLdir/$id.ms$mersize.le$merlimit.mcdat\n";
        if (runCommand($cmd)) {
            unlink "$MERYLdir/$id.ms$mersize.le$merlimit.mcidx";
            unlink "$MERYLdir/$id.ms$mersize.le$merlimit.mcdat";
            die "Failed to count mers lessthanorequal $merlimit in $id\n";
        }
    }

    return "$id.ms$mersize.le$merlimit";
}


#  Return the number of mers in a meryl file.
#
sub numberOfMers ($) {
    my $mers = 0;
    open(F, "$meryl -Dc -s $_[0] |");
    while (<F>) {
        $mers = $1 if (m/Found\s(\d+)\smers/);
    }
    close(F);
    print STDERR "$_[0] has $mers mers.\n";
    return($mers);
}




sub buildMask {

    return if (-e "$ATACdir/$matches.mask.done");

        my $minFile="min.$mercount1.$mercount2";

        #  $mercount1 and $mercount2 are the mers we want to use for
        #  searching.  Obviously, only in-common mers can be found, we
        #  make a file of those mers here.

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

        #  From that list of in-common mers (in-common and below some
        #  count) we want to make a list of the mers that can be used in
        #  the search table.  We can either make a positive (use these
        #  mers) or negative (don't use these mers) list, we just want to
        #  pick the smaller of the two.
        #
        #
        #  The positive 'include' list is just the 'min' mers found above.
        #
        #  The negative 'exclude' list is the min mers, removed from the mers in id1.
        #
        my $includeSize = (-s "$ATACdir/$minFile.mcdat");
        my $excludeSize = (-s "$MERYLdir/$id1.ms$mersize.mcdat") - (-s "$ATACdir/$minFile.mcdat");

        print STDERR "includeSize is proportional to $includeSize.\n";
        print STDERR "excludeSize is proportional to $excludeSize.\n";

        #  But this sometimes breaks (if the mcidx files are different sizes), so we now
        #  pay the cost of actually counting the number of mers.
        #
        $includeSize = numberOfMers("$ATACdir/$minFile");
        $excludeSize = numberOfMers("$MERYLdir/$id1.ms$mersize") - $includeSize;

        print STDERR "includeSize is exactly $includeSize mers.\n";
        print STDERR "excludeSize is exactly $excludeSize mers.\n";

        if ($includeSize < $excludeSize) {
            rename "$ATACdir/$minFile.mcidx", "$ATACdir/$matches.include.mcidx";
            rename "$ATACdir/$minFile.mcdat", "$ATACdir/$matches.include.mcdat";
        } else {
            if (! -e "$ATACdir/$matches.exclude.mcdat") {
                print STDERR "Finding 'exclude' mers!\n";

                #  Our use of xor here is really just a subtraction.  We
                #  want to report those mers that are only in the first
                #  file, not in the second.  All mers in the second file
                #  should be in the first file, by construction.

                my $cmd;
                $cmd  = "$meryl ";
                $cmd .= "-M xor ";
                $cmd .= "-s $MERYLdir/$id1.ms$mersize ";
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

    exit(0) if ($merylOnly == 1);
}



sub findHits {
    my $segmentID   = "000";
    my @segmentIDs;

    open(F, "$leaff -F $MERYLdir/$id1.fasta --partitionmap $numSegments |");
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

    foreach my $segmentID (@segmentIDs) {

        #  For large runs, while developing, we found it very useful
        #  to build the tables first, save them to disk, then do the
        #  compute.  This is also mandatory if one wants to segment
        #  the other assembly to reduce the time each piece runs.
        #
        #  However, doing so adds a lot of complexity to this script,
        #  and isn't terribly useful anymore.
        #

        my $cmd;
        $cmd  = "$seatac ";
        $cmd .= "-verbose ";
        $cmd .= "-mersize     $mersize ";
        $cmd .= "-minlength   $minfill ";
        $cmd .= "-maxgap      $maxgap ";
        $cmd .= "-numthreads  $numThreads ";
        $cmd .= "-table       $MERYLdir/$id1.fasta ";
        $cmd .= "-stream      $MERYLdir/$id2.fasta ";
        $cmd .= "-only        $ATACdir/$matches.include " if (-e "$ATACdir/$matches.include.mcdat");
        $cmd .= "-mask        $ATACdir/$matches.exclude " if (-e "$ATACdir/$matches.exclude.mcdat");
        $cmd .= "-use         $ATACdir/$matches-segment-$segmentID ";
        $cmd .= "-output      $ATACdir/$matches-segment-$segmentID.matches ";
        $cmd .= "-filtername  $filtername " if (defined($filtername));
        $cmd .= "-filteropts  \"-1 $id1 -2 $id2 $filteropts\" ";
        $cmd .= "> $ATACdir/$matches-$segmentID.out 2>&1";

        if (! -e "$ATACdir/$matches-segment-$segmentID.matches") {
            if (runCommand($cmd)) {
                unlink "$ATACdir/$matches-segment-$segmentID.matches.crash";
                rename "$ATACdir/$matches-segment-$segmentID.matches", "$ATACdir/$matches-segment-$segmentID.matches.crash";
                die "Failed to run $matches-$segmentID\n";
            }
        }
    }

    return(@segmentIDs);
}



sub sortMatches (@) {
    my @segmentIDs = @_;

    return if (-e "$ATACdir/$matches.matches.sorted");

    my $mfiles;

    #  Check that each search finished, and build a list of all the match files.
    #
    foreach my $segmentID (@segmentIDs) {
        if (-e "$ATACdir/$matches-segment-$segmentID.matches") {
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

    if ($crossSpecies){
        # The non-default parameters for Mouse versus Rat.
        print ATACFILE "/matchExtenderMinEndRunLen=4\n";
        print ATACFILE "/matchExtenderMaxMMBlock=5\n";
        print ATACFILE "/matchExtenderMinBlockSep=5\n";
        print ATACFILE "/matchExtenderMinIdentity=0.7\n";
        print ATACFILE "/matchExtenderMaxNbrSep=100\n";
        print ATACFILE "/matchExtenderMaxNbrPathMM=25\n";
        print ATACFILE "/globalMatchMinSize=20\n";
        print ATACFILE "/fillIntraRunGapsErate=0.30\n";

        $matchExtenderOpts = "-e 4 -b 5 -s 5 -i 0.70 -p 100 -d 25";
    }

    #print ATACFILE "/matchesFile=$ATACdir/$matches.matches.sorted.bz2\n";

    close(ATACFILE);



    #  The original sort used keys seq1, seq2, pos1, pos2.  The next
    #  consumer of this data, matchextender certainly benefits from
    #  the first two keys, and probably needs the last two....nope, it
    #  does an internal sort given all matches for a pair of
    #  sequences.

    my $tmpdir;
    $tmpdir = "$ATACdir";
    $tmpdir = "/scratch" if (-d "/scratch");

    #  Kill LANG, gnu grep sucks if this is set.
    delete $ENV{'LANG'};

    #  Sort the first index
    #
    #print STDERR "Sorting by the FIRST axis.\n";
    #if (runCommand("cat $mfiles | grep -- \^M\  | sort -y -T $tmpdir -k 5,5 -k 9,9 >> $ATACdir/$matches.matches.sorted")) {
    #    die "Failed to sort $ATACdir!\n";
    #}

    #  Sort the second index, useful if that one is chromosome.
    #
    print STDERR "Sorting by the SECOND axis.\n";
    if (runCommand("cat $mfiles | grep -- '^M ' | sort -y -T $tmpdir -k 9,9 -k 5,5 >> $ATACdir/$matches.matches.sorted")) {
        die "Failed to sort $ATACdir!\n";
    }
}





sub extendMatches {

    return if (-e "$ATACdir/$matches.matches.sorted.extended");

    my $cmd;
    $cmd  = "$BINdir/matchExtender $matchExtenderOpts ";
    $cmd .= "< $ATACdir/$matches.matches.sorted ";
    $cmd .= "> $ATACdir/$matches.matches.sorted.extended ";

    if (runCommand($cmd)) {
        rename("$ATACdir/$matches.matches.sorted.extended",
               "$ATACdir/$matches.matches.sorted.extended.FAILED");
        die "matchExtender failed.\n";
    }
}



sub makeChains {

    return if (-e "$ATACdir/$matches.atac");

    if (!defined($ENV{"TMPDIR"})) {
        print STDERR "WARNING:  TMPDIR not set, defaulting to '$ATACdir'.\n";
        $ENV{"TMPDIR"} = $ATACdir;
    }

    #  Path to the python shared-objects (in lib) and the python scripts.
    #
    $ENV{'PYTHONPATH'} = "$LIBdir";

    if (runCommand("python $chainer $ATACdir/$matches.matches.sorted.extended")) {
        print STDERR "PYTHONPATH=$ENV{'PYTHONPATH'}\n";
        die "Chainer failed.\n";
    }

    if (! -e "$ATACdir/$matches.atac") {
        system("ln -s $ATACdir/$matches.matches.sorted.extended.ckpLast $ATACdir/$matches.atac");
    }
}





#  OK, no, really finish with clumps.
#
sub makeClumps {

    return if (-e "$ATACdir/$matches.atac.clumps5000");

    my $cmd;
    $cmd  = "cd $ATACdir && ";
    $cmd .= "grep ^M $ATACdir/$matches.atac | ";
    $cmd .= "cut -d' ' -f 1-12 | ";
    $cmd .= "sort -k5,5 -k6n ";
    $cmd .= "> $ATACdir/$matches.atac.sorted && ";
    $cmd .= "$BINdir/clumpMaker -c 5000 -2 -S -f $ATACdir/$matches.atac.sorted ";
    $cmd .= "> $ATACdir/$matches.atac.clumps5000";

    if (runCommand($cmd)) {
        rename("$ATACdir/$matches.atac.clumps5000", "$ATACdir/$matches.atac.clumps5000.FAILED");
        die "Failed to make clumps!\n";
    }
}



#  Utility to run a command and check the exit status.  We used to try
#  to decode the exit status...sigh.
#
sub runCommand {
    my $cmd = shift @_;

    print STDERR "$cmd\n";

    if (system($cmd)) {
        return(1);
    }
    return(0);
}
