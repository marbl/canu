#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

$| = 1;

use FindBin;
use lib "$FindBin::Bin";
use libBri;

use strict;

use vars qw($personality $exechome $searchGENOME $mergeCounts $filterMRNA $filterEST $sim4db $leaff $cleanPolishes $toFILTER $sortHits $sortPolishes $parseSNPs $pickBest);

my $scriptVersion = "8";
my $startTime   = time();

$exechome      = "$FindBin::Bin";
$searchGENOME  = "$exechome/searchGENOME";
$mergeCounts   = "$exechome/mergeCounts";
$filterMRNA    = "$exechome/filterMRNA";
$filterEST     = "$exechome/filterEST";
$sim4db        = "$exechome/sim4db";
$leaff         = "$exechome/leaff";
$cleanPolishes = "$exechome/cleanPolishes";
$toFILTER      = "$exechome/filterPolishes";
$sortHits      = "$exechome/sortHits";
$sortPolishes  = "$exechome/sortPolishes";
$parseSNPs     = "$exechome/parseSNP";
$pickBest      = "$exechome/pickBestPolish";
$personality   = "-help";


die "Can't find/execute $searchGENOME\n"  if (! -x $searchGENOME);
die "Can't find/execute $mergeCounts\n"   if (! -x $mergeCounts);
die "Can't find/execute $filterMRNA\n"    if (! -x $filterMRNA);
die "Can't find/execute $filterEST\n"     if (! -x $filterEST);
die "Can't find/execute $sim4db\n"        if (! -x $sim4db);
die "Can't find/execute $leaff\n"         if (! -x $leaff);
die "Can't find/execute $cleanPolishes\n" if (! -x $cleanPolishes);
die "Can't find/execute $sortHits\n"      if (! -x $sortHits);
die "Can't find/execute $sortPolishes\n"  if (! -x $sortPolishes);
die "Can't find/execute $parseSNPs\n"     if (! -x $parseSNPs);
die "Can't find/execute $pickBest\n"      if (! -x $pickBest);


require "util/checkArgs.pl";
require "util/1-configure.pl";
require "util/2-search.pl";
require "util/3-filter.pl";
require "util/4-polish.pl";
require "util/5-assemble.pl";


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
        ($opt eq "-mapmrna")) {
        $personality = $opt;
    }
}


#  If we are asked to restart, ensure that we are doing mapest,
#  mapmrna or mapsnp, read in the original args, and prepend those to
#  the current args.  This lets us override original options, e.g.,
#  the quality levels.
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
        if (($opt eq "-mapest")  ||
            ($opt eq "-mapmrna") ||
            ($opt eq "-mapsnp")) {
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



checkArgs(@ARGV);

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

    $runInformationFile = createRunInformation($dir, "-mapest-nofilter", "$dir", "$est", "$gen", @ARGV);

    configure      ("-configure", "$dir", "-genomic", "$gen", @ARGV);
    search         ("-searchest", "$dir", "-cdna", "$est", "-mersize", "20", "-species", "human", @ARGV);
    filter         ("-filternone", "$dir", @ARGV);
    polish         ("-polish", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
    assembleOutput ("-assembleoutput", "$dir", "-mincoverage", "50", "-minidentity", "95", @ARGV);
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
        system("$parseSNPs $snpdelimiter $snpsizetag $snppostag $snpoffset -F $dir/snps-failed -O $dir/snps-parsed < $dir/polishes-good.sorted-by-cDNA > $dir/summary-snps");
    }
} else {
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

exit(0);
