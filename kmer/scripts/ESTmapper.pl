#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

$| = 1;

use strict;

#  Change Log:
#
#  Tue Apr 30 17:35:29 EDT 2002
#  1) Alignments are reported by default.  Use -noalign to turn them supress them.
#  2) Each search will be attempted four times.  If it fails to complete
#     after four tries, the next search is started.  The entire ESTmapper
#     run will fail after all searches are attempted.
#
#  Wed May  8 14:23:46 EDT 2002
#  1) Added checks to the end of the search phase to ensure that the output
#     is 100% complete.
#  2) Added "-stats" and "-nostats" options to turn on/off run-time statistics.
#     The default is no run-time statistics.
#
#  Thu May  9 10:51:58 EDT 2002
#  Added '-p' flag to all mkdir calls
#
#  Wed May 22 19:45:31 EDT 2002
#  scriptVersion = 2
#  1) Added "-searchthreads" to tell the search to use N threads.  Meaning of
#     "-local" is now the number of concurrent processes to run
#  2) Searches now run concurrently, as specified by "-local".  The search output
#     is checked after all searches have been attempted once.  If any searches
#     have failed, they are restarted, up to three times.
#  3) The default for searches is to run one process with four threads.  For
#     polishing, the default is to run four concurrently.
#  4) -local -> -localpolishes and -localsearches
#     -farm  -> -farmpolishes and -farmsearches
#  5) Times are now reported for polishing
#  6) Added -hitsortmemory to set the memory use when sorting hits.
#     The default is 300MB.
#
#  Mon Jul  1 15:07:18 EDT 2002
#  scriptVersion = 3
#  Added the long intron fix.  It is enabled by default.
#  Added "-nocleanup" and "-cleanup" to disable/enable it.
#
#  Wed Jul 24 15:36:17 EDT 2002
#  scriptVersion = 4
#  1) Added "-longintron" to set the size of a long intron
#  2) Replaced perl match filter with C version.
#
#  Mon Aug 12 13:42:48 EDT 2002
#  Small changes to farm search submission -- it now correctly detects that
#  farm jobs have finished.
#
#  Tue Oct 22 15:54:27 EDT 2002
#  scriptVersion = 5
#  Fixed a low-level problem where non-IUPAC base names would terminate
#  sequences abnormally early when reverse complemented (e.g., the complement
#  of 'X' was \0, which terminated the string at that point).  This
#  affected at least sim4db (rarely).
#
#  Also added DEVELOPMENT support for mapping of SNP's.  Use at your own risk,
#  until it's finalized.
#
#  Tue Dec  3 13:38:25 EST 2002
#  scriptVersion = 6
#
#  1) Added -mapest-nofilter to do no filtering when mapping ests.
#  2) Added -minlength (also to sim4dbseq)
#  3) Added -species to set species specific options (e.g., the mask file)
#  4) Added -maxintron to enable the region splitting on large intron.  Default
#     is 2Mbp.
#  5) Removed polishes-short and polishes-lowquality from the output
#  6) Removed -good and -short, replaced with -minlength, -mincoverage and -minidentity
#


my $scriptVersion = "6";           # analogous to a buildNumber
my $scriptDate    = "03dec02";

#  We use the libBri package of usefull stuff.  It's located in the same place
#  as the ESTmapper.pl script.  That's what FindBin tells us.
#
use FindBin;
use lib "$FindBin::Bin";
use libBri;

my $startTime   = time();

my $exechome      = "$FindBin::Bin";

if ($exechome =~ m!^/cluster(/.*)$!) {
    $exechome = $1;
}
while ($exechome =~ m!^/member\w+(/.*)$!) {
    $exechome = $1;
}

print STDERR "Using exechome of '$exechome'\n";

my $searchGENOME  = "$exechome/searchGENOME";
my $filterMRNA    = "$exechome/filterMRNA";
my $filterEST     = "$exechome/filterEST";
my $sim4db        = "$exechome/sim4db";
my $leaff         = "$exechome/leaff";
my $cleanPolishes = "$exechome/cleanPolishes";
my $toFILTER      = "$exechome/filterPolishes";
my $sortHits      = "$exechome/sortHits";
my $sortPolishes  = "$exechome/sortPolishes";
my $parseSNPs     = "$exechome/parseSNP";

die "Can't find/execute $searchGENOME\n"  if (! -x $searchGENOME);
die "Can't find/execute $filterMRNA\n"    if (! -x $filterMRNA);
die "Can't find/execute $filterEST\n"     if (! -x $filterEST);
die "Can't find/execute $sim4db\n"        if (! -x $sim4db);
die "Can't find/execute $leaff\n"         if (! -x $leaff);
die "Can't find/execute $cleanPolishes\n" if (! -x $cleanPolishes);
die "Can't find/execute $sortHits\n"      if (! -x $sortHits);
die "Can't find/execute $sortPolishes\n"  if (! -x $sortPolishes);
die "Can't find/execute $parseSNPs\n"     if (! -x $parseSNPs);

sub showHelp {
    print STDERR "Basic help:\n";
    print STDERR "\n";
    print STDERR "  ESTmapper.pl -mapest  work-directory ests.fasta genomic.fasta\n";
    print STDERR "\n";
    print STDERR "or\n";
    print STDERR "\n";
    print STDERR "  ESTmapper.pl -mapmrna work-directory mrna.fasta genomic.fasta\n";
    print STDERR "\n";
    print STDERR "Read the manual for details.\n";
}



#  Create a .runInformaiton file, containing supposedly useful information
#  about this run.
#
sub createRunInformation {
    my ($dir, @ARGS) = @_;
    my $time = time();

    system("mkdir -p $dir") if (! -d "$dir");

    my $runInformationFile = "$dir/.runInformation.$time";

    print STDERR "ESTmapper: This is run $runInformationFile\n";
    print STDERR "ESTmapper: The following options are in effect:";
    foreach my $opt (@ARGS) {
        print STDERR "\nESTmapper:\t" if ($opt =~ m/^-/);
        print STDERR "$opt ";
    }
    print STDERR "\n";

    #  Write some information and the args to a run info file
    #
    open(F, "> $runInformationFile");
    print F "startTime:  $time (", scalar(localtime($time)), ")\n";
    print F "operator:   $ENV{'USER'}\n";
    print F "host:       " . `uname -a`;
    print F "version:    $scriptVersion\n";
    print F "parameters:";
    foreach my $a (@ARGS) { print F " $a"; }
    print F "\n";
    close(F);

    unlink "$dir/.runInformation";
    system("ln -s $dir/.runInformation.$time $dir/.runInformation");

    #  Write the current set of args to the runOptions file
    #
    open(F, "> $dir/.runOptions");
    foreach my $a (@ARGS) { print F "$a\n"; }
    close(F);

    return($runInformationFile);
}




my $personality = "-help";
foreach my $opt (@ARGV) {
    if (($opt eq "-restart")         ||
        ($opt eq "-configure")       ||
        ($opt eq "-searchest")       ||
        ($opt eq "-searchmrna")      ||
        ($opt eq "-filterest")       ||
        ($opt eq "-filtermrna")      ||
        ($opt eq "-filternone")      ||
        ($opt eq "-polish")          ||
        ($opt eq "-assembleoutput")  ||
        ($opt eq "-mapest")          ||
        ($opt eq "-mapest-nofilter") ||
        ($opt eq "-mapsnp")          ||
        ($opt eq "-mapesttoest")     ||
        ($opt eq "-mapmrna")) {
        $personality = $opt;
    }
}


#
#  Restart only supported for -mapest and -mapmrna.
#
if ($personality eq "-restart") {
    shift @ARGV;
    my $dir = shift @ARGV;

    if (! -e "$dir/.runOptions") {
        print STDERR "ESTmapper/restart-- Nothing to restart!\n";
        exit;
    }

    open(F, "< $dir/.runOptions");
    my @origARGS = <F>;
    chomp @origARGS;
    close(F);

    @ARGV = (@origARGS, @ARGV);

    $personality = "-help";
    foreach my $opt (@ARGV) {
        if (($opt eq "-mapest") || ($opt eq "-mapmrna")) {
            $personality = $opt;
        }
    }

    print STDERR "ESTmapper: Reastarting with the following options:";
    foreach my $opt (@ARGV) {
        print STDERR "\nESTmapper:\t" if ($opt =~ m/^-/);
        print STDERR "$opt ";
    }
    print STDERR "\n";
}


my $runInformationFile;


if ($personality eq "-mapest") {
    shift @ARGV;
    my $dir = shift @ARGV;
    my $est = shift @ARGV;
    my $gen = shift @ARGV;

    $runInformationFile = createRunInformation($dir, "-mapest", "$dir", "$est", "$gen", @ARGV);

    configure      ("-configure", "$dir", "-genomic", "$gen", @ARGV);
    search         ("-searchest", "$dir", "-cdna", "$est", "-mersize", "20", "-species", "human", @ARGV);
    filter         ("-filterest", "$dir", @ARGV);
    polish         ("-polish", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
    assembleOutput ("-assembleoutput", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
} elsif ($personality eq "-mapest-nofilter") {
    shift @ARGV;
    my $dir = shift @ARGV;
    my $est = shift @ARGV;
    my $gen = shift @ARGV;

    $runInformationFile = createRunInformation($dir, "-mapest", "$dir", "$est", "$gen", @ARGV);

    configure      ("-configure", "$dir", "-genomic", "$gen", @ARGV);
    search         ("-searchest", "$dir", "-cdna", "$est", "-mersize", "20", "-species", "human", @ARGV);
    filter         ("-filternone", "$dir", @ARGV);
    polish         ("-polish", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
    assembleOutput ("-assembleoutput", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
} elsif ($personality eq "-mapsnp") {
    shift @ARGV;
    my $dir = shift @ARGV;
    my $est = shift @ARGV;
    my $gen = shift @ARGV;

    $runInformationFile = createRunInformation($dir, "-mapsnp", "$dir", "$est", "$gen", @ARGV);

    configure      ("-configure", "$dir", "-genomic", "$gen", @ARGV);
    search         ("-searchsnp", "$dir", "-cdna", "$est", "-mersize", "20", "-species", "human", @ARGV);
    filter         ("-filtersnp", "$dir", @ARGV);
    polish         ("-polish", "$dir", "-mincoverage", "80", "-minidentity", "95", @ARGV);
    assembleOutput ("-assembleoutput", "$dir", "-mincoverage", "80", "-minidentity", "95", @ARGV);

    #  Parse the SNPs out
    #
    if (! -e "$dir/snps-parsed") {
        print STDERR "ESTmapper--  Parsing the SNPs\n";

        #  Sort, if needed.
        #
        if (! -e "$dir/polishes-good.sorted-by-cDNA") {
            print STDERR "ESTmapper--  Sorting polishes by sequence ID; using 4GB memory maximum.\n";
            system("$sortPolishes -m 4000 -c -v < $dir/polishes-good > $dir/polishes-good.sorted-by-cDNA");
        }

        #  Parse the options, looking for SNP specific ones
        #
        my @ARGS = @ARGV;
        my $snpdelimiter = "";
        my $snpsizetag   = "";
        my $snppostag    = "";
        my $snpoffset    = "";

        while (scalar @ARGS > 0) {
            my $arg = shift @ARGS;

            if ($arg eq "-snpdelimiter") {
                $arg = shift @ARGS;
                $snpdelimiter = "-d \"$arg\"";
            } elsif ($arg eq "-snpsizetag") {
                $arg = shift @ARGS;
                $snpsizetag = "-s \"$arg\"";
            } elsif ($arg eq "-snppostag") {
                $arg = shift @ARGS;
                $snppostag = "-p \"$arg\"";
            } elsif ($arg eq "-snpoffset") {
                $arg = shift @ARGS;
                $snpoffset = "-o $arg";
            }
        }

        #  PARSE!
        #
        print ("$parseSNPs $snpdelimiter $snpsizetag $snppostag $snpoffset -F $dir/snps-failed -O $dir/snps-parsed < $dir/polishes-good.sorted-by-cDNA\n");
        system("$parseSNPs $snpdelimiter $snpsizetag $snppostag $snpoffset -F $dir/snps-failed -O $dir/snps-parsed < $dir/polishes-good.sorted-by-cDNA");
    }

} elsif ($personality eq "-mapmrna") {
    shift @ARGV;
    my $dir = shift @ARGV;
    my $est = shift @ARGV;
    my $gen = shift @ARGV;

    $runInformationFile = createRunInformation($dir, "-mapmrna", "$dir", "$est", "$gen", @ARGV);

    configure      ("-configure", "$dir", "-genomic", "$gen", @ARGV);
    search         ("-searchmrna", "$dir", "-cdna", "$est", "-mersize", "20", "-species", "human", @ARGV);
    filter         ("-filtermrna", "$dir", @ARGV);
    polish         ("-polish", "$dir", "-mincoverage", "50", "-minidentity", "95", "-relink", "1000", "-abort", @ARGV);
    assembleOutput ("-assembleoutput", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
} elsif ($personality eq "-mapesttoest") {
    shift @ARGV;
    my $dir = shift @ARGV;
    my $est = shift @ARGV;
    my $gen = shift @ARGV;

    configure      ("-configure", "$dir", "-genomic", "$gen", @ARGV);
    search         ("-searchest", "$dir", "-cdna", "$est", "-mersize", "20", "-species", "human", @ARGV);
    filter         ("-filterest", "$dir", @ARGV);

    #  Fix the hits to remove all the me-to-me hits
    #
    if (! -e "$dir/all/2-filter/filteredHits.self") {
        print STDERR "Removing self hits.\n";
        rename "$dir/all/2-filter/filteredHits", "$dir/all/2-filter/filteredHits.self";
        open(F, "< $dir/all/2-filter/filteredHits.self");
        open(G, "> $dir/all/2-filter/filteredHits");
        while (<F>) {
            if (m/-.\s-e\s(\d+)\s-D\s(\d+)/) {
                print G $_ if ($1 ne $2);
            } else {
                print STDERR "Hit error!\n";
            }
        }
        close(G);
        close(F);
    }

    polish         ("-polish", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
    assembleOutput ("-assembleoutput", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
} else {
    #createRunInformation($dir, @ARGV);

    showHelp()               if ($personality eq "-help");
    configure(@ARGV)         if ($personality eq "-configure");
    search(@ARGV)            if ($personality eq "-searchest");
    search(@ARGV)            if ($personality eq "-searchmrna");
    search(@ARGV)            if ($personality eq "-searchsnp");
    filter(@ARGV)            if ($personality eq "-filtermrna");
    filter(@ARGV)            if ($personality eq "-filterest");
    polish(@ARGV)            if ($personality eq "-polish");
    assembleOutput(@ARGV)    if ($personality eq "-assembleoutput");
}

print STDERR "ESTmapper: script finished everything in ", time() - $startTime, " wall-clock seconds.\n" if (time() != $startTime);


if (defined($runInformationFile)) {
    my $time = time();

    open(F, ">> $runInformationFile");
    print F "endTime:  $time (", scalar(localtime($time)), ")\n";
    close(F);
}

exit;


#  Make sure that all the polishes are finished and OK.
#  Returns 0 if all done, 1 if not done.
#
sub polishesNotDone {
    my ($path) = @_;

    my $failed = 0;

    open(F, "< $path/3-polish/run-script");
    while (!eof(F)) {
        my $idx = <F>;  chomp $idx;
        my $cmd  = <F>;
        if (! -e "$path/3-polish/$idx.touch") {
            print STDERR "Polish $idx failed to finish.\n";
            $failed = 1;
        }
    }
    close(F);

    return $failed;
}


sub assembleOutput {
    my $startTime  = time();
    my $errHdr     = "ERROR: ESTmapper/assemble--";
    my @ARGS       = @_;
    my $path       = "";
    my $minc       = 50;
    my $mini       = 95;
    my $minl       = 0;
    my $doCleaning =  1;
    my $longthresh =  100000;

    print STDERR "ESTmapper: Performing an assembleOutput.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-assembleoutput") {
            $path = shift @ARGS;
        }

        if ($arg eq "-mincoverage") {
            $minc = shift @ARGS;
        }
        if ($arg eq "-minidentity") {
            $mini = shift @ARGS;
        }
        if ($arg eq "-minlength") {
            $minl = shift @ARGS;
        }

        if ($arg eq "-cleanup") {
            $doCleaning = 1;
        }
        if ($arg eq "-nocleanup") {
            $doCleaning = 0;
        }
        if ($arg eq "-longintron") {
            $doCleaning = 1;
            $longthresh = int(shift @ARGS)
        }
    }

    ($path eq "") and die "$errHdr no directory given.\n";
    (! -d "$path") and die "$errHdr no directory '$path' found!\n";
    (($mini < 0) || ($mini > 100)) and die "$errHdr supply a value 0 <= x <= 100 for minidentity!\n";
    (($minc < 0) || ($minc > 100)) and die "$errHdr supply a value 0 <= x <= 100 for mincoverage!\n";
    ($minl < 0) and die "$errHdr supply a value x >= 0 for minlength!\n";

    (polishesNotDone($path) > 0) and die "There are unfinished polishing jobs.\n";

    #  Check that the filtering is compatable with the polishing.
    #
    if (-e "$path/3-polish/parameters") {
        open(F, "< $path/3-polish/parameters");
        $_ = <F>;
        $_ = <F>;

        my $miniL = int(<F>);  #  Quality values used for last filtering
        my $mincL = int(<F>);
        my $minlL = int(<F>);

        my $miniP = int(<F>);  #  Quality values used for polishing
        my $mincP = int(<F>);
        my $minlP = int(<F>);
        close(F);

        if ($mini < $miniP) {
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Percent identity quality level too low for existing polishing!\n";
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Polished at percent align-sequence identity = %3d, requested filtration at %3d.\n", $miniP, $mini;
        }
        if ($minc < $mincP) {
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Coverage quality level too low for existing polishing!\n";
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Polished at percent query-sequence identity = %3d, requested filtration at %3d.\n", $mincP, $minc;
        }
        if ($minl < $minlP) {
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Length quality level too low for existing polishing!\n";
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Polished at length = %3d, requested filtration at %3d.\n", $minlP, $minl;
        }

        #  If the filter quality has changed, we need to refilter.  Nuke
        #  the filterLevel file, print a message.
        #
        if (($mini != $miniL) ||
            ($minc != $mincL) ||
            ($minl != $minlL)) {
            print STDERR "ESTmapper/assembleOutput-- filtering criteria changed; refiltering.\n";

            printf STDERR "ESTmapper/assembleOutput-- identity:  percent align-sequence identity: old=%3d new=%3d\n", $miniL, $mini;
            printf STDERR "ESTmapper/assembleOutput-- coverage:  percent query-sequence identity: old=%3d new=%3d\n", $mincL, $minc;
            printf STDERR "ESTmapper/assembleOutput-- length:    length in bp of match:           old=%3d new=%3d\n", $minlL, $minl;

            unlink "$path/polishes-good";
            unlink "$path/polishes-goodshort";
            unlink "$path/polishes-lowquality";
            unlink "$path/summary";
        }
    } else {
        die "ESTmapper/assemblyOutput-- ERROR: Couldn't find polishing parameters.  Script error.\n";
    }



    if (! -e "$path/polishes-good") {
        print STDERR "ESTmapper/assembleOutput-- filtering polishes by quality.\n";

        print STDERR "ESTmapper/assembleOutput-- identity:  percent align-sequence identity: $mini\n";
        print STDERR "ESTmapper/assembleOutput-- coverage:  percent query-sequence identity: $minc\n";
        print STDERR "ESTmapper/assembleOutput-- length:    length in bp of match:           $minl\n";

        #  Take all the polishes and filter them into
        #  polishes-good, polishes-goodshort, polishes-lowquality
        #
        my $cmd;

        #$cmd  = "find $path/3-polish/ -name '*.polished' -print | sort | xargs -n 100 cat | ";
        #$cmd .= "$cleanPolishes -threshold $longthresh -savejunk | " if ($doCleaning);
        #$cmd .= "$toFILTER -c $shrt -i $good -o $path/polishes-good -j $path/polishes-aborted | ";
        #$cmd .= "$toFILTER -c 0     -i $good -o $path/polishes-goodshort > $path/polishes-lowquality";

        $cmd  = "find $path/3-polish/ -name '*.polished' -print | sort | xargs -n 100 cat | ";
        $cmd .= "$cleanPolishes -threshold $longthresh -savejunk | " if ($doCleaning);
        $cmd .= "$toFILTER -c $minc -i $mini -l $minl -o $path/polishes-good -j $path/polishes-aborted > /dev/null";

        system($cmd);

        unlink "$path/cDNA-good.fasta";
        unlink "$path/cDNA-missing.fasta";
        unlink "$path/cDNA-repeat.fasta";
        unlink "$path/cDNA-zero.fasta";
        unlink "$path/summary";
    } else {
        print STDERR "ESTmapper/assembleOutput-- polishes already filtered.\n";
    }

    #
    #  Segregate the sequences
    #

    if (! -e "$path/cDNA-missing.fasta") {
        unlink "$path/summary";

        #  For each polished file, figure out the ESTs that belong with those polishes.
        #  cDNA-good.fasta, cDNA-goodshort.fasta, cDNA-lowquality.fasta, cDNA-zero.fasta, cDNA-missing.fasta
        #
        print STDERR "ESTmapper/assembleOutput-- finding 'good' cDNA.\n";
        &libBri::splitFastABasedOnPolishes("$path/0-input/cDNA.fasta",
                                           "$path/polishes-good",
                                           "$path/cDNA-good.fasta",
                                           "$path/cDNA-lost.fasta");

        #print STDERR "ESTmapper/assembleOutput-- finding 'good' cDNA.\n";
        #&libBri::splitFastABasedOnPolishes("$path/0-input/cDNA.fasta",
        #                                   "$path/polishes-good",
        #                                   "$path/cDNA-good.fasta",
        #                                   "$path/cDNA-good.fasta.negative");

        #print STDERR "ESTmapper/assembleOutput-- finding 'good, but short' cDNA.\n";
        #&libBri::splitFastABasedOnPolishes("$path/cDNA-good.fasta.negative",
        #                                   "$path/polishes-goodshort",
        #                                   "$path/cDNA-goodshort.fasta",
        #                                   "$path/cDNA-goodshort.fasta.negative");

        #print STDERR "ESTmapper/assembleOutput-- finding 'low quality' cDNA.\n";
        #&libBri::splitFastABasedOnPolishes("$path/cDNA-goodshort.fasta.negative",
        #                                   "$path/polishes-lowquality",
        #                                   "$path/cDNA-lowquality.fasta",
        #                                   "$path/cDNA-lost.fasta");

        if (-e "$path/2-filter/repeats") {
            print STDERR "ESTmapper/assembleOutput-- finding 'repeat' cDNA.\n";
            system("$leaff -I $path/0-input/cDNA.fasta -q $path/2-filter/repeats > $path/cDNA-repeat.fasta");
        }

        print STDERR "ESTmapper/assembleOutput-- finding 'zero hit' cDNA.\n";
        &libBri::copyZeroFastA("$path/0-input/cDNA.fasta",
                               "$path/2-filter/hitCounts",
                               "$path/cDNA-zero.fasta");

        #  subtractFastAfromFastA checks for input file existance, so we need
        #  not worry that cDNA-repeat.fasta might not exist
        #
        print STDERR "ESTmapper/assembleOutput-- finding 'missing' cDNA.\n";
        &libBri::subtractFastAfromFastA("$path/cDNA-lost.fasta",
                                        "$path/cDNA-repeat.fasta",
                                        "$path/cDNA-zero.fasta",
                                        "$path/cDNA-missing.fasta");

        #  Remove the temporary files
        #
        #unlink "$path/cDNA-good.fasta.negative";
        #unlink "$path/cDNA-goodshort.fasta.negative";
        unlink "$path/cDNA-lost.fasta";
    }

    #
    #  Summarize
    #

    if ((! -e "$path/summary") || (-z "$path/summary")) {
        my ($mat, $est, $scf);

        open(F, "> $path/summary");

        #  Summarize the polished sets
        #
        print STDERR "ESTmapper/assembleOutput-- counting 'good' matches.\n";
        ($mat, $est, $scf) = &libBri::summarizePolishes("$path/polishes-good");
        print F "GOOD: >= $mini% identity, >= $minc% composite, >= $minl bp\n"; 
        if ($mat > 0) {
            print F "cDNA-genomic matches  $mat matches ($est different cDNA and $scf genomic)\n";
            print F "Matches per cDNA      ", int(10000 * $mat / $est) / 10000.0, " matches/cDNA\n";
            print F "Matches per genomic   ", int(10000 * $mat / $scf) / 10000.0, " matches/genomic\n";
        } else {
            print F "cDNA-genomic matches  None.\n";
        }
        print F "\n";

        #print STDERR "ESTmapper/assembleOutput-- counting 'good, but short' matches.\n";
        #($mat, $est, $scf) = &libBri::summarizePolishes("$path/polishes-goodshort");
        #print F "GOOD but SHORT: >= $good% identity, < $shrt% composite\n"; 
        #if ($mat > 0) {
        #    print F "cDNA-genomic matches  $mat matches ($est different cDNA and $scf genomic)\n";
        #    print F "Matches per cDNA      ", int(10000 * $mat / $est) / 10000.0, " matches/cDNA\n";
        #    print F "Matches per genomic   ", int(10000 * $mat / $scf) / 10000.0, " matches/genomic\n";
        #} else {
        #    print F "cDNA-genomic matches  None.\n";
        #}
        #print F "\n";

        #print STDERR "ESTmapper/assembleOutput-- counting 'all the good' matches.\n";
        #my ($mat, $est, $scf) = &libBri::summarizePolishes("$path/polishes-good", "$path/polishes-goodshort");
        #print F "ALL THE GOOD:  (both 'GOOD' and 'GOOD but SHORT')\n";
        #if ($mat > 0) {
        #    print F "cDNA-genomic matches  $mat matches ($est different cDNA and $scf genomic)\n";
        #    print F "Matches per cDNA      ", int(10000 * $mat / $est) / 10000.0, " matches/cDNA\n";
        #    print F "Matches per genomic   ", int(10000 * $mat / $scf) / 10000.0, " matches/genomic\n";
        #} else {
        #    print F "cDNA-genomic matches  None.\n";
        #}
        #print F "\n";

        #print STDERR "ESTmapper/assembleOutput-- counting 'low quality' matches.\n";
        #($mat, $est, $scf) = &libBri::summarizePolishes("$path/polishes-lowquality");
        #print F "LOW-QUALITY: < $good% identity\n";
        #if ($mat > 0) {
        #    print F "cDNA-genomic matches  $mat matches ($est different cDNA and $scf genomic)\n";
        #    print F "Matches per cDNA      ", int(10000 * $mat / $est) / 10000.0, " matches/cDNA\n";
        #    print F "Matches per genomic   ", int(10000 * $mat / $scf) / 10000.0, " matches/genomic\n";
        #} else {
        #    print F "cDNA-genomic matches  None.\n";
        #}
        #print F "\n";

        print STDERR "ESTmapper/assembleOutput-- counting cDNA.\n";
        print F "cDNA COUNTS:\n";
        my $cnttotl = int(`grep -c '>' $path/0-input/cDNA.fasta`);
        my $cntgood = int(`grep -c '>' $path/cDNA-good.fasta`);
        #my $cntshrt = int(`grep -c '>' $path/cDNA-goodshort.fasta`);
        #my $cntlowq = int(`grep -c '>' $path/cDNA-lowquality.fasta`);
        my $cntmiss = int(`grep -c '>' $path/cDNA-missing.fasta`);
        my $cntrept = int(`grep -c '>' $path/cDNA-repeat.fasta`) if (-e "$path/cDNA-repeat.fasta");
        my $cntzero = int(`grep -c '>' $path/cDNA-zero.fasta`);

        printf F "cDNA:            %8d\n", $cnttotl, "\n";
        printf F "cDNA-good:       %8d (%8.4f%%)\n", $cntgood, 100 * $cntgood / $cnttotl;
        #printf F "cDNA-goodshort:  %8d (%8.4f%%)\n", $cntshrt, 100 * $cntshrt / $cnttotl;
        #printf F "cDNA-lowquality: %8d (%8.4f%%)\n", $cntlowq, 100 * $cntlowq / $cnttotl;
        printf F "cDNA-missing:    %8d (%8.4f%%)\n", $cntmiss, 100 * $cntmiss / $cnttotl;
        printf F "cDNA-repeat:     %8d (%8.4f%%)\n", $cntrept, 100 * $cntrept / $cnttotl  if (-e "$path/cDNA-repeat.fasta");
        printf F "cDNA-zero:       %8d (%8.4f%%)\n", $cntzero, 100 * $cntzero / $cnttotl;
    }

    print STDERR "ESTmapper: assembleOutput script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}



sub polish {
    my $startTime = time();
    my $errHdr     = "ERROR: ESTmapper/polish--";
    my @ARGS       = @_;
    my $path       = "";
    my $mini       = "95";
    my $minc       = "50";
    my $minl       = "0";

    my $minsim4i   = "90";
    my $minsim4c   = "45";
    my $minsim4l   = "0";

    my $always     = "";
    my $relink     = "";
    my $batchsize  = 0;
    my $numbatches = 256;
    my $farm       = 0;
    my $farmqueue  = "";
    my $farmcode   = "";
    my $local      = 1;
    my $numcpus    = 4;
    my $runnow     = 1;
    my $aligns     = "-align";
    my $stats      = 0;
    my $abort      = "";

    print STDERR "ESTmapper: Performing a polish.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-polish") {
            $path = shift @ARGS;
        }

        if ($arg eq "-mincoverage") {
            $minc = shift @ARGS;
        }
        if ($arg eq "-minidentity") {
            $mini = shift @ARGS;
        }
        if ($arg eq "-minlength") {
            $minl = shift @ARGS;
        }

        if ($arg eq "-minsim4coverage") {
            $minsim4c = shift @ARGS;
        }
        if ($arg eq "-minsim4identity") {
            $minsim4i = shift @ARGS;
        }
        if ($arg eq "-minsim4length") {
            $minsim4l = shift @ARGS;
        }

        if ($arg eq "-alwaysprint") {
            $always = "-alwaysprint " . shift @ARGS;
        }
        if ($arg eq "-relink") {
            $relink = "-H " . shift @ARGS;
        }
        if ($arg eq "-batchsize") {
            $batchsize = int(shift @ARGS);
        }
        if ($arg eq "-numbatches") {
            $numbatches = int(shift @ARGS);
        }
        if ($arg eq "-farmpolishes") {
            $farm     = 1;
            $local    = 0;
            $farmqueue = shift @ARGS;
            $farmcode  = shift @ARGS
        }
        if ($arg eq "-localpolishes") {
            $farm  = 0;
            $local = 1;
            $numcpus = shift @ARGS;
        }
        if ($arg eq "-runlater") {
            $runnow = 0;
        }
        if ($arg eq "-aligns") {
            $aligns = "-align";
        }
        if ($arg eq "-noaligns") {
            $aligns = "";
        }
        if ($arg eq "-stats") {
            $stats = 1;
        }
        if ($arg eq "-nostats") {
            $stats = 0;
        }
        if ($arg eq "-abort") {
            $abort = "-Mp 0.25 -Ma 10000";
        }
    }

    ($path eq "") and die "$errHdr no directory given.\n";
    (! -d "$path") and die "$errHdr no directory '$path' found!\n";

    system("mkdir $path/3-polish") if (! -d "$path/3-polish");


    #  Save the parameters, these are used on later invocations of
    #  polish, and in filter to make sure the user isn't an idiot.
    #
    if (-e "$path/3-polish/parameters") {
        print STDERR "ESTmapper/polish-- Using original parameters.\n";

        open(F, "< $path/3-polish/parameters");
        $numbatches = int(<F>);
        $batchsize  = int(<F>);
        $mini       = <F>;      chomp $mini;
        $minc       = <F>;      chomp $minc;
        $minl       = <F>;      chomp $minl;
        $minsim4i   = <F>;      chomp $minsim4i;
        $minsim4c   = <F>;      chomp $minsim4c;
        $minsim4l   = <F>;      chomp $minsim4l;
        $relink     = <F>;      chomp $relink;
        $always     = <F>;      chomp $always;
        $aligns     = <F>;      chomp $aligns;
        $abort      = <F>;      chomp $abort;
        close(F);

        print STDERR "ESTmapper/polish-- Polish quality suitable for $minsim4i percent identity and\n";
        print STDERR "ESTmapper/polish--                             $minsim4c percent coverage\n";
        print STDERR "ESTmapper/polish-- To rerun polishes at a different quality level,\n";
        print STDERR "ESTmapper/polish-- remove the 3-polish directory.\n";
    } else {

        #  Do a little error checking; if both $batchsize and
        #  $numbatches are zero, set $batchsize to make 256 batches.
        #
        if (($batchsize == 0) && ($numbatches == 0)) {
            $numbatches = 256;
        }

        #  If $batchsize is not specified, compute it.
        #
        if ($batchsize == 0) {
            $batchsize = int(`wc -l < $path/2-filter/filteredHits` / $numbatches) + 1;
            $batchsize = 500 if ($batchsize < 500);
        }

        #  Adjust the sim4 qualities based on the final quality desired
        #
        $mini = 0 if ($mini < 0);
        $minc = 0 if ($minc < 0);
        $minl = 0 if ($minl < 0);

        $minsim4i = $mini - 5 if ($mini - 5 < $minsim4i);
        $minsim4c = $minc - 5 if ($minc - 5 < $minsim4c);
        $minsim4l = $minl     if ($minl     < $minsim4l);

        $minsim4i = 0 if ($minsim4i < 0);
        $minsim4c = 0 if ($minsim4c < 0);
        $minsim4l = 0 if ($minsim4l < 0);

        #  Save the parameters
        #
        open(F, "> $path/3-polish/parameters");
        print F "$numbatches\n$batchsize\n";
        print F "$mini\n$minc\n$minl\n";
        print F "$minsim4i\n$minsim4c\n$minsim4l\n";
        print F "$relink\n$always\n$aligns\n$abort\n";
        close(F);
    }


    #  Display what parameters we are using
    #
    print STDERR "ESTmapper/polish--   minidentity = $mini ($minsim4i)\n";
    print STDERR "ESTmapper/polish--   mincoverage = $minc ($minsim4c)\n";
    print STDERR "ESTmapper/polish--   minlength   = $minl ($minsim4l)\n";
    print STDERR "ESTmapper/polish--   relink      = $relink\n";
    print STDERR "ESTmapper/polish--   always      = $always\n";
    print STDERR "ESTmapper/polish--   aligns      = $aligns\n";
    print STDERR "ESTmapper/polish--   abort       = $abort\n";


    #  Splits the filteredHits into several pieces, and outputs a script
    #  that runs sim4db on those pieces.
    #
    if (! -e "$path/3-polish/run-script") {
        print STDERR "ESTmapper/polish-- Creating scripts with $batchsize lines in each.\n";

        my $idx = "0000";

        open(H, "< $path/2-filter/filteredHits");
        open(S, "> $path/3-polish/run-script");
        while (!eof(H)) {
            my $c = 0;

            open(F, "> $path/3-polish/$idx.scr");
            while (($c < $batchsize) && (!eof(H))) {
                $_ = <H>;
                print F $_;
                $c++;
            }
            close(F);

            #  The run-script is composed of three lines per command:
            #    the first line is the index of the run
            #    the second is the command
            #
            print S "$idx\n";
            print S "bsub -q $farmqueue -o $path/3-polish/$idx.stdout -R \"select[physmem>300]rusage[physmem=600]\" -P $farmcode " if ($farm);
            print S "$sim4db -cdna $path/0-input/cDNA.fasta -genomic $path/0-input/genomic.fasta ";
            print S "$aligns $always $relink $abort -cut 0.6 ";
            print S "-mincoverage $minsim4c ";
            print S "-minidentity $minsim4i ";
            print S "-minlength $minsim4l ";
            print S "-script $path/3-polish/$idx.scr ";
            print S "-output $path/3-polish/$idx.polished ";
            print S "-stats  $path/3-polish/$idx.stats " if ($stats == 1);
            print S "-touch  $path/3-polish/$idx.touch\n";

            $idx++;
        }
        close(S);
        close(H);

        print STDERR "ESTmapper/polish-- Created $idx scripts.\n";
    }

    #  Builds a list of things to run by looking at the run-script, and
    #  seeing if the output for a given command exists.
    #
    my $polishesToPerform = 0;

    open(F, "< $path/3-polish/run-script");
    open(S, "> $path/3-polish/run.sh");
    while (!eof(F)) {
        my $idx = <F>;  chomp $idx;
        my $cmd = <F>;

        if (! -e "$path/3-polish/$idx.touch") {
            $polishesToPerform = 1;
            print S $cmd;
        }
    }
    close(S);
    close(F);

    #  Wipe any summaries, cDNA-* and polished files if we need to polish more stuff.
    #
    if ($polishesToPerform == 1) {
        print STDERR "ESTmapper/polish-- more polishes to compute - removing old output.\n";
        unlink "$path/cDNA-good.fasta";
        unlink "$path/cDNA-goodshort.fasta";
        unlink "$path/cDNA-lowquality.fasta";
        unlink "$path/cDNA-missing.fasta";
        unlink "$path/cDNA-repeat.fasta";
        unlink "$path/cDNA-zero.fasta";
        unlink "$path/polishes-aborted";
        unlink "$path/polishes-good";
        unlink "$path/polishes-goodshort";
        unlink "$path/polishes-lowquality";
        unlink "$path/summary";
    }


    #  Run things, or tell the user to do it for us.
    #
    if ($runnow) {
        if ($local && $polishesToPerform) {
            print STDERR "ESTmapper/polish-- Running locally, $numcpus at a time.\n";

            &libBri::schedulerSetNumberOfProcesses($numcpus);
            &libBri::schedulerSetShowCommands(0);
            &libBri::schedulerSetShowStatus(1);
            open(F, "< $path/3-polish/run.sh");
            while (<F>) {
                chomp;
                &libBri::schedulerSubmit($_);
            }
            close(F);
            &libBri::schedulerFinish();
        }

        if ($farm && $polishesToPerform) {
            print STDERR "ESTmapper/polish-- Submitting to the farm.\n";

            open(F, "< $path/3-polish/run.sh");
            while (<F>) {
                chomp;
                system($_);
            }
            close(F);
            
            #  Hang around waiting for them to finish.
            die "ESTmapper/polish-- I don't know how to monitor LSF jobs!\nESTmapper/polish-- Please restart when they're all done.\n";
        }

        unlink "$path/3-polish/run.sh";
    } else {
        die "ESTmapper/polish-- Please run the stuff in\nESTmapper/polish--   $path/3-polish/run.sh\n";
    }

    my $clkTime = 0;
    my $sysTime = 0;
    my $usrTime = 0;
    open(F, "< $path/3-polish/run-script");
    while (!eof(F)) {
        my $idx = <F>;  chomp $idx;
        my $cmd = <F>;

        if (-e "$path/3-polish/$idx.stats") {
            open(X, "< $path/3-polish/$idx.stats");
            while (<X>) {
                if (m/^clockTime:\s+(\d+\.\d+)$/) {
                    $clkTime += $1;
                }
                if (m/^systemTime:\s+(\d+\.\d+)$/) {
                    $sysTime += $1;
                }
                if (m/^userTime:\s+(\d+\.\d+)$/) {
                    $usrTime += $1;
                }
            }
            close(X);
        }
    }
    close(F);

    print STDERR "ESTmapper: sim4db required $clkTime seconds wall-clock time.\n";
    print STDERR "ESTmapper: sim4db required $sysTime seconds system time.\n";
    print STDERR "ESTmapper: sim4db required $usrTime seconds user time.\n";
    print STDERR "ESTmapper: Polish script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


sub readScafList {
    my $f = shift @_;
    open(F, "< $f");
    my @l = <F>;
    close(F);
    chomp @l;
    return(@l);
}

sub filter {
    my $startTime = time();
    my $errHdr    = "ERROR: ESTmapper/filter--";
    my @ARGS      = @_;
    my $path      = "";
    my $type      = "";
    my $hitMemory = "600";
    my $verbose   = "";

    print STDERR "ESTmapper: Performing a filter.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-filterest") {
            $path = shift @ARGS;
            $type  = "est";
        }
        if ($arg eq "-filtersnp") {
            $path = shift @ARGS;
            $type  = "snp";
        }
        if ($arg eq "-filtermrna") {
            $path = shift @ARGS;
            $type  = "mrna";
        }
        if ($arg eq "-filternone") {
            $path = shift @ARGS;
            $type  = "none";
        }
        if ($arg eq "-hitsortmemory") {
            $hitMemory = shift @ARGS;
        }
        if ($arg eq "-verbose") {
            $verbose = "-verbose";
        }
    }

    ($path eq "") and die "$errHdr no directory given.\n";
    (! -d "$path") and die "$errHdr no directory '$path' found!\n";
    (! -f "$path/0-input/scaffolds-list") and die "$errHdr no scaffolds-list?\n";
    (! -f "$path/1-search/allDone") and die "The searches are not complete.\n";

    system("mkdir $path/2-filter") if (! -d "$path/2-filter");

    my @scafList = readScafList("$path/0-input/scaffolds-list");

    #
    #  Merge all the hit counts into one list
    #
    if (! -e "$path/2-filter/hitCounts") {
        my @counts;

        print STDERR "ESTmapper/search-- Merging counts.\n";
        foreach my $s (@scafList) {
            my $idx = 0;

            open(F, "< $path/1-search/$s.count");
            while (<F>) {
                $counts[$idx++] += int($_);
            }
            close(F);
        }

        print STDERR "ESTmapper/search-- Writing counts.\n";
        open(F, "> $path/2-filter/hitCounts");
        foreach my $c (@counts) {
            print F "$c\n";
        }
        close(F);
    }

    if (! -e "$path/2-filter/filteredHits") {
        my $fcmd;
        my $scmd;

        if ($type eq "est") {
            #print STDERR "\n\nWARNING:  USING MODIFIED VALUE FOR REPEAT THRESHOLD.  SHOULD BY 100!\n\n\n";
            $fcmd  = "$filterEST $verbose -u 100 -r 100 -q 0.2 -c $path/2-filter/hitCounts ";
            $fcmd .= "-Fu $path/2-filter/uniqHits ";
            $fcmd .= "-Ff $path/2-filter/filtHits ";
            $fcmd .= "-Fr $path/2-filter/repeats ";
            $fcmd .= "$path/1-search/*hits";

            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter/sorttmp $path/2-filter/uniqHits $path/2-filter/filtHits > $path/2-filter/filteredHits";
        }

        if ($type eq "snp") {
            $fcmd  = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits > $path/2-filter/filtHits";

            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter/sorttmp $path/2-filter/filtHits > $path/2-filter/filteredHits";
        }

        if ($type eq "mrna") {
            $fcmd  = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits > $path/2-filter/filtHits";

            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter/sorttmp $path/2-filter/filtHits > $path/2-filter/filteredHits";
        }

        if ($type eq "none") {
            $fcmd  = "";

            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter/sorttmp $path/1-search/*hits > $path/2-filter/filteredHits";

            print STDERR "RUNNING FILTER: $scmd\n";
        }

        print STDERR "ESTmapper/filter-- Filtering.\n";
        system($fcmd);

        print STDERR "ESTmapper/filter-- Sorting.\n";
        system($scmd);
   }

    print STDERR "ESTmapper: Filter script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}





sub countTheCountFile {
    my ($path, $s) = @_;
    my $c = 0;

    #  If the hits file is NOT found, remove the count file.
    #
    if (! -e "$path/1-search/$s.hits") {
        unlink "$path/1-search/$s.count";
    }

    #  If there is a count file, count the number of lines in it.
    #  There should be one for each cDNA.
    #
    if (-e "$path/1-search/$s.count") {
        $c = int(`wc -l < $path/1-search/$s.count`);
    }

    print STDERR "ESTmapper/search-- found $c cDNA for search $s\n";

    return($c);
}



sub search {
    my $startTime = time();
    my $errHdr     = "ERROR: ESTmapper/search--";
    my @ARGS       = @_;
    my $path       = "";
    my $opts       = "";
    my $cdna       = "";
    my $mersize    = 20;
    my $maskFile   = "";
    my $verbose    = "";
    my $stats      = 0;
    my $farm       = 0;
    my $farmqueue  = "";
    my $farmcode   = "";
    my $local      = 1;
    my $numthread  = 2;
    my $numproc    = 4;
    my $runnow     = 1;
    my $maxintron  = "-maxintron 2000000";
    my $species    = "";

    print STDERR "ESTmapper: Performing a search.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-searchest") {
            $path = shift @ARGS;
            $opts = "-singlelength 20 -multiplelength 30 -smallsequence 100";
        }
        if ($arg eq "-searchmrna") {
            $path = shift @ARGS;
            $opts = "-singlelength 30 -multiplelength 50 -smallsequence 0";
        }
        if ($arg eq "-searchsnp") {
            $path = shift @ARGS;
            $opts = "-singlecoverage 0.3 -multiplecoverage 0.3 -smallsequence 10000000 -extendminimum 100 -extendweight 2";
        }
        if ($arg eq "-searchopts") {
            $opts .= " " . shift @ARGS;
            print "$opts\n";
        }
        if ($arg eq "-cdna") {
            $cdna = shift @ARGS;
        }
        if ($arg eq "-mersize") {
            $mersize = shift @ARGS;
        }
        if ($arg eq "-maskmers") {
            $species  = "";
            $maskFile = shift @ARGS;
        }
        if ($arg eq "-nomaskmers") {
            $maskFile = "";
        }
        if ($arg eq "-verbose") {
            $verbose = "-verbose";
        }
        if ($arg eq "-stats") {
            $stats = 1;
        }
        if ($arg eq "-nostats") {
            $stats = 0;
        }
        if ($arg eq "-farmsearches") {
            $farm     = 1;
            $local    = 0;
            $farmqueue = shift @ARGS;
            $farmcode  = shift @ARGS
        }
        if ($arg eq "-localsearches") {
            $farm  = 0;
            $local = 1;
            $numproc = shift @ARGS;
        }
        if ($arg eq "-searchthreads") {
            $numthread = shift @ARGS;
        }
        if ($arg eq "-runlater") {
            $runnow = 0;
        }
        if ($arg eq "-maxintron") {
            $maxintron = "-maxintron " . shift @ARGS;
        }
        if ($arg eq "-species") {
            $species  = shift @ARGS;
            $maskFile = "";
            if ($species eq "none") {
                $maskFile = "";
            } elsif (($species eq "human") ||
                     ($species eq "hum")) {
                $maskFile = "$exechome/data/frequentMers-hum-20.fasta";
            } elsif (($species eq "mouse") ||
                     ($species eq "mus")) {
                $maskFile = "$exechome/data/frequentMers-mus-20.fasta";
            } elsif (($species eq "rattus") ||
                     ($species eq "rat")) {
                $maskFile = "$exechome/data/frequentMers-rat-20.fasta";
            } else {
                print STDERR "ESTmapper/search-- Unknown species '$species'.\n";
                $maskFile = "";
            }
        }
    }

    ($path eq "") and die "$errHdr no directory given.\n";
    (! -d "$path") and die "$errHdr no directory '$path' found!\n";
    (! -f "$path/0-input/scaffolds-list") and die "$errHdr no scaffolds-list?\n";
    (($maskFile ne "") && (! -e "$maskFile")) and die "$errHdr Can't find -maskmers '$maskFile'!\n";
    ($cdna eq "") and die "$errHdr Please supply some cDNA to search for!\n";
    (! -f $cdna) and die "$errHdr cDNA sequences '$cdna' not found!\n";


    if ($maskFile eq "") {
        print STDERR "ESTmapper/search-- No file of maskmers - mer masking disabled.\n";
    } else {
        print STDERR "ESTmapper/search-- Using '$maskFile' for mer masking.\n";
    }

    system("mkdir $path/1-search") if (! -d "$path/1-search");

    if (-e "$path/1-search/allDone") {
        print STDERR "ESTmapper/search-- this stage all done, and results checked previously.\n";
        return;
    }

    #  Save a pointer to the cDNA used
    #
    if (! -e "$path/0-input/cDNA.fasta") {
        system("ln -s $cdna      $path/0-input/cDNA.fasta");
    }

    if (! -e "$path/0-input/cDNA.fastaidx") {
        if (-f "${cdna}idx") {
            system("ln -s ${cdna}idx $path/0-input/cDNA.fastaidx");
        } else {
            print STDERR "Building index for cDNA sequences.\n";
            system("$leaff -I $path/0-input/cDNA.fasta");
        }
    }

    my @scafList = readScafList("$path/0-input/scaffolds-list");

    #  Create a bunch of scripts to process
    #
    #  Rewrite the command everytime.  This fixes the problem where
    #  we would, say, change the stats or number of threads.
    #
    foreach my $s (@scafList) {
        open(F, "> $path/1-search/$s.cmd");
        print F "$searchGENOME $verbose -binary -mersize $mersize $opts $maxintron -numthreads $numthread";
        print F " -cdna $path/0-input/cDNA.fasta";
        print F " -genomic $path/0-input/genomic.fasta";
        print F " -use $path/0-input/scaffolds-$s";
        print F " -output $path/1-search/$s.hits";
        print F " -count $path/1-search/$s.count";
        print F " -stats $path/1-search/$s.stats" if ($stats == 1);
        print F " -mask $maskFile" if ($maskFile ne "");
        close(F);
        system("chmod 755 $path/1-search/$s.cmd");
    }

    my $cdnaInInput  = int(`$leaff -F $path/0-input/cDNA.fasta -d`);

    #  Run searches.  If the search terminated properly, the
    #  hit-counts file should exist.  Run (maybe re-run) the search if
    #  it isn't there.
    #
    if ($runnow) {
        if ($local) {
            print STDERR "ESTmapper/search-- Local mode requested; ", scalar @scafList, " processes to compute,\n";
            print STDERR "ESTmapper/search-- Local mode requested; $numproc concurrent processes,\n";
            print STDERR "ESTmapper/search-- Local mode requested; each with $numthread threads.\n";

            #  Build a list of things to run.
            #
            my @searchesToRun;

            foreach my $s (@scafList) {
                if (countTheCountFile($path, $s) != $cdnaInInput) {
                    push @searchesToRun, $s;
                }
            }

            my $numTries = 0;
          again:
            $numTries++;

            #  Run the searches.  We use the scheduler, then check
            #  everything at the end.  This is a little less friendly
            #  to the user, but much easier for the implementor.
            #
            if (scalar(@searchesToRun) > 0) {
                &libBri::schedulerSetNumberOfProcesses($numproc);
                &libBri::schedulerSetShowCommands(0);
                &libBri::schedulerSetShowStatus(1);
                foreach my $s (@searchesToRun) {
                    &libBri::schedulerSubmit("sh $path/1-search/$s.cmd");
                }
                &libBri::schedulerFinish();

                #  See if anything failed.
                #
                undef @searchesToRun;
                print STDERR "ESTmapper/search-- checking search output.  All should have $cdnaInInput cDNA.\n";
                foreach my $s (@scafList) {
                    if (countTheCountFile($path, $s) != $cdnaInInput) {
                        push @searchesToRun, $s;
                    }
                }
            }

            if (($numTries < 3) && (scalar(@searchesToRun) > 0)) {
                print STDERR "ESTmapper/search-- ", scalar(@searchesToRun), " searches failed.  Retrying.\n";
                goto again;
            }

            if (scalar(@searchesToRun) > 0) {
                print STDERR "ESTmapper/search-- Searches failed.\n";
                foreach my $s (@searchesToRun) {
                    print STDERR "ESTmapper/search-- Search $s failed.  Output saved as *.CRASH\n";
                    rename "$path/1-search/$s.count", "$path/1-search/$s.count.CRASH";
                    rename "$path/1-search/$s.hits", "$path/1-search/$s.hits.CRASH";
                }
            }
        }

        if ($farm) {
            open(F, "< $path/0-input/memoryLimit");
            my $farmMemory = <F>; chomp $farmMemory;
            close(F);

            print STDERR "ESTmapper/search-- Farm mode requested; ", scalar @scafList, " processes to compute,\n";
            print STDERR "ESTmapper/search-- Farm mode requested; each with $numthread threads,\n";
            print STDERR "ESTmapper/search-- Farm mode requested; $farmMemory MB per process.\n";

            for (my $x=10; $x>0; $x--) {
                print STDERR "ESTmapper/search-- Abort in $x seconds to change. \r";
                sleep(1);
            }
            print STDERR "ESTmapper/search-- STARTING FARM RUN!           \n";

            my $jobsToRun = 0;

            foreach my $s (@scafList) {
                if (countTheCountFile($path, $s) != $cdnaInInput) {
                    print STDERR "ESTmapper/search-- submitting search $s\n";
                    #print STDERR "bsub -q $farmqueue -o $path/1-search/$s.stdout -R \"select[physmem>$farmMemory]rusage[physmem=$farmMemory]\" -P $farmcode $path/1-search/$s.cmd";
                    system("bsub -q $farmqueue -o $path/1-search/$s.stdout -R \"select[physmem>$farmMemory]rusage[physmem=$farmMemory]\" -P $farmcode $path/1-search/$s.cmd");
                    $jobsToRun++;
                } else {
                    print STDERR "ESTmapper/search-- search $s finished successfully!\n";
                }
            }

            #  Hang around waiting for them to finish.
            if ($jobsToRun > 0) {
                die "ESTmapper/search-- I don't know how to monitor LSF jobs!\nESTmapper/search-- Please restart when they're all done.\n";
            }
        }
    }


    #  See if anything is still broken.  If it did, holler!
    #  We also sum the time used here.
    #
    my $died    = 0;
    my $sysTime = 0;
    my $usrTime = 0;
    foreach my $s (@scafList) {
        if (! -e "$path/1-search/$s.count") {
            print STDERR "$errHdr SEARCH $s DIED!\n";
            $died++;
        }
        if (-e "$path/1-search/$s.stats") {
            open(F, "< $path/1-search/$s.stats");
            while (<F>) {
                if (m/^systemTime:\s+(\d+\.\d+)$/) {
                    $sysTime += $1;
                }
                if (m/^userTime:\s+(\d+\.\d+)$/) {
                    $usrTime += $1;
                }
            }
            close(F);
        }
    }
    die "Searches have a problem ($died).\n" if ($died > 0);


    #  Rather lazy way to inform the next step (and future calls to this step) that we're all done.
    #
    system("touch $path/1-search/allDone");

    print STDERR "ESTmapper: searchGENOME required $sysTime seconds system time.\n";
    print STDERR "ESTmapper: searchGENOME required $usrTime seconds user time.\n";
    print STDERR "ESTmapper: Search script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}



sub configure {
    my $startTime = time();
    my $errHdr    = "ERROR: ESTmapper/configure--";
    my @ARGS      = @_;
    my $path      = "";
    my $genomic   = "";
    my $memory    = 800;
    my $seqlen    = 0;

    print STDERR "ESTmapper: Performing a configure.\n";

    #
    #  Parse the args to find the path, then read any saved
    #  configuration, then reparse the args.
    #

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-configure") {
            $path = shift @ARGS;
        }
        if ($arg eq "-genomic") {
            $genomic = shift @ARGS;
        }
        if ($arg eq "-memory") {
            $memory = shift @ARGS;
        }
    }

    #  400000000 bases -> 3500M -->  8.75 (old)
    #  125828998 bases -> size=2400M res=1500M --> 11.92 (from 03dec02)
    #  104857600 bases -> 
    #
    my $scaleFactor = 12;

    #  Determine the length of each piece
    #
    $seqlen = $memory * 1024.0 * 1024.0 / $scaleFactor;

    ($path eq "") and die "$errHdr no directory given.\n";
    ($genomic eq "") and die "$errHdr no genomic sequence given.\n";
    (! -f $genomic) and die "$errHdr can't find the genomic sequence '$genomic'\n";

    #  Make a place for us to work
    #
    system("mkdir $path")          if (! -d "$path");
    system("mkdir $path/0-input")  if (! -d "$path/0-input");
    system("mkdir $path/1-search") if (! -d "$path/1-search");
    system("mkdir $path/2-filter") if (! -d "$path/2-filter");
    system("mkdir $path/3-polish") if (! -d "$path/3-polish");

    #  For the farm, we need to save the amount of memory we can use.
    #
    open(F, "> $path/0-input/memoryLimit") or die "Can't write $path/0-input/memoryLimit\n";
    print F "$memory\n";
    close(F);

    if (! -e "$path/0-input/scaffolds-list") {
        print STDERR "ESTmapper/configure-- Use about ${memory}MB -> $seqlen bases per chunk.\n";

        #  Remember the genomic file for later
        #
        system("ln -s ${genomic}    $path/0-input/genomic.fasta")    if ((-e "${genomic}")    && (! -e "$path/0-input/genomic.fasta"));
        system("ln -s ${genomic}idx $path/0-input/genomic.fastaidx") if ((-e "${genomic}idx") && (! -e "$path/0-input/genomic.fastaidx"));
        system("ln -s ${genomic}inf $path/0-input/genomic.fastainf") if ((-e "${genomic}inf") && (! -e "$path/0-input/genomic.fastainf"));

        if (! -f "$path/0-input/genomic.fasta") {
            die "$errHdr can't find the genomic sequence '$path/0-input/genomic.fasta'\n";
        }

        if ((! -f "$path/0-input/genomic.fastaidx") && (! -f "$path/0-input/genomic.fastainf")) {
            print STDERR "ESTmapper/configure-- Generating the index and info for '$path/0-input/genomic.fasta'\n";
            print STDERR "ESTmapper/configure-- WARNING:  This is done in the work directory!\n";
            system("$leaff -F $path/0-input/genomic.fasta -i > $path/0-input/genomic.fastainf");
        }

        if (! -f "$path/0-input/genomic.fastaidx") {
            print STDERR "ESTmapper/configure-- Generating the index for '$path/0-input/genomic.fasta'\n";
            print STDERR "ESTmapper/configure-- WARNING:  This is done in the work directory!\n";
            system("$leaff -I $path/0-input/genomic.fasta");
        }
        if (! -f "$path/0-input/genomic.fastainf") {
            print STDERR "ESTmapper/configure-- Generating the info for '$path/0-input/genomic.fasta'\n";
            print STDERR "ESTmapper/configure-- WARNING:  This is done in the work directory!\n";
            system("$leaff -F $path/0-input/genomic.fasta -i > $path/0-input/genomic.fastainf");
        }

        #
        #  Partition the genome into itty-bitty pieces
        #

        my @scaffolds;

        #  The original method (for original leaff)
        #
        #open(F, "< $path/0-input/genomic.fastainf");
        #while (<F>) {
        #    if (m/^\s*(\d+)\]\s+\d+\s+(\d+)$/) {
        #        push @scaffolds, "$2.$1";
        #    } else {
        #        die "$errHdr Invalid line in information file: $_.\n";
        #    }
        #}
        #close(F);

        my $idx = 0;
        open(F, "< $path/0-input/genomic.fastainf");
        while (<F>) {
            #  Skip all the fasta file information
            last if (m/^\d+\s+\d+\s+(\d+)\s+\d+\D/);
        }
        while (<F>) {
            if (m/^\d+\s+\d+\s+(\d+)\s+\d+\D/) {
                push @scaffolds, "$1.$idx";
                $idx++;
            } else {
                die "$errHdr Invalid line in information file: $_.\n";
            }
        }
        close(F);

        @scaffolds = reverse sort { $a <=> $b } @scaffolds;

        open(L, "> $path/0-input/scaffolds-list");

        my $outputFile = "000";

        #  Special case; for any scaffolds more than the limit, write
        #  them now.
        #
        my $l;
        my $i;

        ($l, $i) = split '\.', $scaffolds[0];

        while ((scalar @scaffolds > 0) &&
               ($l >= $seqlen)) {

            print L "$outputFile\n";

            open(F, "> $path/0-input/scaffolds-$outputFile");
            print F "$i\n";
            close(F);

            printf STDERR "ESTmapper/configure-- WARNING:  Scaffold $i ($l bases) requires %8.3fMB of memory!\n", $l * $scaleFactor / 1024.0 / 1024.0;

            shift @scaffolds;

            ($l, $i) = split '\.', $scaffolds[0];

            $outputFile++;
        }


        #  Now, pack the "small" scaffolds together.
        #
        while (scalar @scaffolds > 0) {
            my $totalLength = 0;
            my @leftover;
            
            undef @leftover;

            print L "$outputFile\n";

            open(F, "> $path/0-input/scaffolds-$outputFile");
            foreach my $v (@scaffolds) {
                my ($l, $i) = split '\.', $v;

                if ($totalLength + $l < $seqlen) {
                    print F "$i\n";
                    $totalLength += $l;
                } else {
                    push @leftover, $v;
                }
            }
            close(F);

            printf STDERR "ESTmapper/configure-- Created group with $totalLength bases (%8.3fMB of memory).\n", $totalLength * $scaleFactor / 1024.0 / 1024.0;

            $outputFile++;
            
            @scaffolds = @leftover;
        }

        close(L);
    }

    print STDERR "ESTmapper: Configure script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}
