
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

package canu::SequenceStore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getNumberOfReadsInStore
             getNumberOfBasesInStore
             getExpectedCoverage
             getSequenceStoreStats
             getSizeOfSequenceStore
             createSequenceStore
             generateReadLengthHistogram);

use strict;
use warnings "all";
no  warnings "uninitialized";

use Cwd qw(getcwd);

use canu::Defaults;
use canu::Execution;

use canu::Report;
use canu::Output;

use canu::Grid_Cloud;


use Carp qw(longmess cluck confess);


#  Return the number of reads in the store.
#
#  The 'all' category returns the number of reads that are in the store; essentailly, the maxID of any read.
#  The 'cor' category returns the number of reads that are available for correction.
#  The 'obt' category returns the number of reads that are available for trimming.
#  The 'utg' category returns the number of reads that are available for assembly.
#
sub getNumberOfReadsInStore ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $nr     = 0;

    confess  if (($tag ne "all") &&
                 ($tag ne "hap") &&
                 ($tag ne "cor") &&
                 ($tag ne "obt") &&
                 ($tag ne "utg"));

    return($nr)   if (! -e "./$asm.seqStore/info.txt");   #  No file, no reads.

    open(F, "< ./$asm.seqStore/info.txt") or caExit("can't open './$asm.seqStore/info.txt' for reading: $!", undef);
    while (<F>) {
        if (m/^\s*(\d+)\s+([0123456789-]+)\s+(.*)\s*$/) {
            $nr = $1 if (($3 eq "total-reads")       && ($tag eq "all"));
            $nr = $1 if (($3 eq "raw")               && ($tag eq "cor"));
            $nr = $1 if (($3 eq "raw")               && ($tag eq "hap"));
            $nr = $1 if (($3 eq "corrected")         && ($tag eq "obt"));
            $nr = $1 if (($3 eq "corrected-trimmed") && ($tag eq "utg"));
        }
    }
    close(F);

    return($nr);
}



#  Return the number of bases in the store.
#
#  Two modes of operation, depending on if the second option is supplied.
#   - If $tag is "cor", "obt" or "utg", return the number of bases in reads
#     that will be used for that stage.
#   - If $tag is not supplied, return the maximum number of bases in
#     any category ('raw', 'raw-trimmed', etc).
#
sub getNumberOfBasesInStore ($@) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $nb     = 0;

    confess  if ((defined($tag)) &&
                 ($tag ne "cor") &&
                 ($tag ne "obt") &&
                 ($tag ne "utg"));

    return($nb)   if (! -e "./$asm.seqStore/info.txt");   #  No file, no bases.

    open(F, "< ./$asm.seqStore/info.txt") or caExit("can't open './$asm.seqStore/info.txt' for reading: $!", undef);
    while (<F>) {
        if (m/^\s*(\d+)\s+(\d+)\s+(.*)\s*$/) {
            $nb = $2   if (($3 eq "raw")               && ($tag eq "cor"));
            $nb = $2   if (($3 eq "corrected")         && ($tag eq "obt"));
            $nb = $2   if (($3 eq "corrected-trimmed") && ($tag eq "utg"));

            $nb = $2   if (($3 eq "raw")               && ($nb < $2) && (!defined($tag)));
            $nb = $2   if (($3 eq "raw-trimmed")       && ($nb < $2) && (!defined($tag)));
            $nb = $2   if (($3 eq "corrected")         && ($nb < $2) && (!defined($tag)));
            $nb = $2   if (($3 eq "corrected-trimmed") && ($nb < $2) && (!defined($tag)));
        }
    }
    close(F);

    return($nb);
}



#  Returns the estimated coverage in reads for the supplied tag (see
#  getNumberOfBasesInStore for an important detail).
#
#  Returns a number rounded to hundredths of a genome, ###.##.
#
sub getExpectedCoverage ($@) {
    my $asm    = shift @_;
    my $tag    = shift @_;

    confess  if ((defined($tag)) &&
                 ($tag ne "cor") &&
                 ($tag ne "obt") &&
                 ($tag ne "utg"));

    return(int(100 * getNumberOfBasesInStore($asm, $tag) / getGlobal("genomeSize")) / 100);
}



sub getSequenceStoreStats ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();

    my ($numRaw, $numRawTri, $numCor, $numCorTri) = (0, 0, 0, 0);
    my ($numPacBio, $numNanopore, $numHiFi)       = (0, 0, 0);

    #  Recreate our metadata dumps, if needed.

    if (! -d "./$asm.seqStore") {
        return(0, 0, 0, 0, 0, 0, 0);
    }

    if (! -e "./$asm.seqStore/info.txt") {
        if (runCommandSilently(".", "$bin/sqStoreDumpMetaData -S ./$asm.seqStore -stats > ./$asm.seqStore/info.txt 2> /dev/null", 1)) {
            caExit("failed to generate $asm.seqStore/info.txt", undef);
        }
    }

    if (! -e "./$asm.seqStore/libraries.txt") {
        if (runCommandSilently(".", "$bin/sqStoreDumpMetaData -S ./$asm.seqStore -libs > ./$asm.seqStore/libraries.txt 2> /dev/null", 1)) {
            caExit("failed to generate $asm.seqStore/libraries.txt", undef);
        }
    }

    #  Count the number of reads or each type.

    open(L, "< ./$asm.seqStore/info.txt") or caExit("can't open './$asm.seqStore/info.txt' for reading: $!", undef);
    while (<L>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        $numRaw    += $v[1]   if (($v[2] eq "raw")               && ($v[1] > 0));
        $numRawTri += $v[1]   if (($v[2] eq "raw-trimmed")       && ($v[1] > 0));
        $numCor    += $v[1]   if (($v[2] eq "corrected")         && ($v[1] > 0));
        $numCorTri += $v[1]   if (($v[2] eq "corrected-trimmed") && ($v[1] > 0));
    }
    close(L);

    open(L, "< ./$asm.seqStore/libraries.txt") or caExit("can't open './$asm.seqStore/libraries.txt' for reading: $!", undef);
    while (<L>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        $numPacBio++       if ($v[1] eq "PacBio");
        $numNanopore++     if ($v[1] eq "Nanopore");
        $numHiFi++         if ($v[1] eq "PacBioHiFi");
    }
    close(L);

    return($numRaw, $numRawTri, $numCor, $numCorTri, $numPacBio, $numNanopore, $numHiFi);
}



sub getSizeOfSequenceStore ($) {
    my $asm    = shift @_;
    my $size   = 0;
    my $idx    = "0000";

    $size += -s "./$asm.seqStore/info";
    $size += -s "./$asm.seqStore/libraries";
    $size += -s "./$asm.seqStore/reads";

    while (-e "./$asm.seqStore/blobs.$idx") {
        $size += -s "./$asm.seqStore/blobs.$idx";
        $idx++;
    }

    return(int($size / 1024 / 1024 / 1024 + 1.5));
}




sub createSequenceStoreMetaDataFiles ($) {
    my $asm = shift @_;
    my $bin = getBinDirectory();

    if (! -e "./$asm.seqStore/info.txt") {
        if (runCommandSilently(".", "$bin/sqStoreDumpMetaData -S ./$asm.seqStore -stats > ./$asm.seqStore/info.txt 2> /dev/null", 1)) {
            caExit("failed to generate $asm.seqStore/info.txt", undef);
        }
    }

    if (! -e "./$asm.seqStore/libraries.txt") {
        if (runCommandSilently(".", "$bin/sqStoreDumpMetaData -S ./$asm.seqStore -libs > ./$asm.seqStore/libraries.txt 2> /dev/null", 1)) {
            caExit("failed to generate $asm.seqStore/libraries.txt", undef);
        }
    }
}


sub createSequenceStore ($@) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    #  If a seqStore exists, make sure that the metadata dumps exist, then
    #  return silently.

    if (-e "./$asm.seqStore") {
        createSequenceStoreMetaDataFiles($asm);
        return;
    }

    #  Fail if there are no inputs.

    caExit("no input files specified, and store not already created, I have nothing to work on!", undef)
        if (scalar(@inputs) == 0);

    #  If any of the files are links to objects, fetch the object
    #  to local disk and update the name.
    #
    #  A similar blcok is used in SequenceStore.pm and HaplotypeReads.pm (twice).

    for (my $ff=0; $ff < scalar(@inputs); $ff++) {
        my ($inType, $inPath) = split '\0', $inputs[$ff];

        if ($inPath =~ m/dnanexus:(.*)=(.*)/) {
            my $link = $1;
            my $name = $2;

            print STDERR "-- Fetch input file './$name' from object '$link'\n";

            fetchFileFromLink($link, $name);

            $inputs[$ff] = "$inType\0./$name";
        }
    }

    if (! -e "./$asm.seqStore.sh") {
        open(F, "> ./$asm.seqStore.sh") or caExit("cant' open './$asm.seqStore.sh' for writing: $0", undef);

        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode(".");
        print F "\n";

        print F "\n";
        print F "$bin/sqStoreCreate \\\n";
        print F "  -o ./$asm.seqStore.BUILDING \\\n";
        print F "  -minlength "  . getGlobal("minReadLength")        . " \\\n";

        if (getGlobal("maxInputCoverage") > 0) {
            print F "  -genomesize " . getGlobal("genomeSize")       . " \\\n";
            print F "  -coverage   " . getGlobal("maxInputCoverage") . " \\\n";
            print F "  -bias       " . getGlobal("readSamplingBias") . " \\\n";
        }

        if (getGlobal("homoPolyCompress") > 0) {   #  Tell sqStoreCreate that these reads should
            print F "  -homopolycompress \\\n";    #  be homopoly compressed, and so length filtering
        }                                          #  should be based on that length.

        foreach my $iii (@inputs) {
            if ($iii =~ m/^(.*)\0(.*)$/) {
                my $tech = $1;
                my $file = $2;
                my @name = split '/', $2;
                my $name = $name[scalar(@name)-1];

                $name = $1   if ($name =~ m/(.*).[xgb][z]2{0,1}$/i);
                $name = $1   if ($name =~ m/(.*).fast[aq]$/i);
                $name = $1   if ($name =~ m/(.*).f[aq]$/i);

                print F "  $tech $name $file \\\n";

            } else {
                caExit("unrecognized sqStoreCreate input file '$iii'", undef);
            }
        }

        print F "&& \\\n";
        print F "mv ./$asm.seqStore.BUILDING ./$asm.seqStore \\\n";
        print F "&& \\\n";
        print F "exit 0\n";
        print F "\n";
        print F "exit 1\n";

        close(F);

        makeExecutable("./$asm.seqStore.sh");
        stashFile("./$asm.seqStore.sh");
    }

    #  Load the store.

    if (runCommand(".", "./$asm.seqStore.sh > ./$asm.seqStore.err 2>&1") > 0) {
        caExit("sqStoreCreate failed; boom!", "./$asm.seqStore.err");
    }

    if (! -e "./$asm.seqStore/info") {
        caExit("sqStoreCreate failed; no 'info' file", "./$asm.seqStore.err");
    }

    createSequenceStoreMetaDataFiles($asm);

    #  Set or unset the homopolymer compression flag.  This is also initially
    #  created by sqStoreCreate (also based on the homoPolyCompress flag) but
    #  we allow canu to explicitly change later.  Though doing that isn't
    #  recommended.

    if (defined(getGlobal("homoPolyCompress"))) {
        touch("./$asm.seqStore/homopolymerCompression");
    } else {
        unlink("./$asm.seqStore/homopolymerCompression");
    }

    #  We don't really know what read type was loaded, and so we attempt to
    #  generate all read length histograms.  If there is no data, no
    #  histogram will be made.

    generateReadLengthHistogram("cor", $asm);
    generateReadLengthHistogram("obt", $asm);
    generateReadLengthHistogram("utg", $asm);

    #  Check if coverage is too low.  We don't actually know (easily) what
    #  type of reads were just loaded, and so we check the maximum found
    #  in the store.

    my $curcov = getExpectedCoverage($asm);
    my $mincov = getGlobal("minInputCoverage");

    if ($curcov < $mincov) {
        print STDERR "--\n";
        print STDERR "-- ERROR:  Read coverage ($curcov) lower than allowed.\n";
        print STDERR "-- ERROR:    minInputCoverage  = $mincov\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  This could be caused by an incorrect genomeSize.\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  You can force Canu to continue by decreasing parameter\n";
        print STDERR "-- ERROR:  minInputCoverage.  Be warned that the quality of corrected\n";
        print STDERR "-- ERROR:  reads and/or contiguity of contigs will be poor.\n";
        print STDERR "--\n";

        caExit("", undef);
    }

    #  Stop if requested.

    stopAfter("sequenceStore");
}



sub generateReadLengthHistogram ($$) {
    my $tag      = shift @_;
    my $asm      = shift @_;
    my $bin      = getBinDirectory();

    my $reads    = getNumberOfReadsInStore($asm, $tag);
    my $bases    = getNumberOfBasesInStore($asm, $tag);
    my $coverage = getExpectedCoverage($asm, $tag);
    my $hist;

    return   if ($reads == 0);

    my @rl;
    my @hi;

    #  Generate a lovely PNG histogram.

    if (! -e "./$asm.seqStore/readlengths-$tag.dat") {
        my $cmd;

        $cmd  = "$bin/sqStoreDumpMetaData \\\n";
        $cmd .= "  -S ./$asm.seqStore \\\n";
        $cmd .= "  -raw \\\n"                  if ($tag eq "cor");
        $cmd .= "  -corrected \\\n"            if ($tag eq "obt");
        $cmd .= "  -corrected -trimmed \\\n"   if ($tag eq "utg");
        $cmd .= "  -uncompressed \\\n";
        $cmd .= "  -lengths \\\n";
        $cmd .= "> ./$asm.seqStore/readlengths-$tag.dat \\\n";
        $cmd .= "2> ./$asm.seqStore/readlengths-$tag.err \n";

        if (runCommandSilently(".", $cmd, 1)) {
            caExit("sqStoreDumpMetaData failed", "./$asm.seqStore/readlengths-$tag.err");
        }

        unlink "./$asm.seqStore/readlengths-$tag.err";

        my $gnuplot = getGlobal("gnuplot");
        my $format  = getGlobal("gnuplotImageFormat");

        open(F, "> ./$asm.seqStore/readlengths-$tag.gp") or caExit("can't open './$asm.seqStore/readlengths-$tag.gp' for writing: $!", undef);
        print F "set title 'read length'\n";
        print F "set xlabel 'read length, bin width = 250'\n";
        print F "set ylabel 'number of reads'\n";
        print F "\n";
        print F "binwidth=250\n";
        print F "set boxwidth binwidth\n";
        print F "bin(x,width) = width*floor(x/width) + binwidth/2.0\n";
        print F "\n";
        print F "set terminal $format size 1024,1024\n";
        print F "set output './$asm.seqStore/readlengths-$tag.$format'\n";
        print F "plot [] './$asm.seqStore/readlengths-$tag.dat' using (bin(\$1,binwidth)):(1.0) smooth freq with boxes title ''\n";
        close(F);

        #  Run gnuplot, but don't care if it blows up.
        if ($gnuplot) {
            runCommandSilently(".", "$gnuplot < /dev/null ./$asm.seqStore/readlengths-$tag.gp > /dev/null 2>&1", 0);
        }
    }

    #  Generate a lovely ASCII histogram.

    if (! -e "./$asm.seqStore/readlengths-$tag.txt") {
        my $cmd;

        $cmd  = "$bin/sqStoreDumpMetaData \\\n";
        $cmd .= "  -S ./$asm.seqStore \\\n";
        $cmd .= "  -raw \\\n"                  if ($tag eq "cor");
        $cmd .= "  -corrected \\\n"            if ($tag eq "obt");
        $cmd .= "  -corrected -trimmed \\\n"   if ($tag eq "utg");
        $cmd .= "  -uncompressed \\\n";
        $cmd .= "  -histogram \\\n";
        $cmd .= "> ./$asm.seqStore/readlengths-$tag.txt \\\n";
        $cmd .= "2> ./$asm.seqStore/readlengths-$tag.err \n";

        if (runCommandSilently(".", $cmd, 1)) {
            caExit("sqStoreDumpMetaData failed", "./$asm.seqStore/readlengths-$tag.err");
        }

        unlink "./$asm.seqStore/readlengths-$tag.err";
    }

    #  Read the ASCII histogram, append to report.

    $hist  = "--\n";
    $hist .= "-- In sequence store './$asm.seqStore':\n";
    $hist .= "--   Found $reads reads.\n";
    $hist .= "--   Found $bases bases ($coverage times coverage).\n";

    if (-e "./$asm.seqStore/readlengths-$tag.txt") {
        open(F, "< ./$asm.seqStore/readlengths-$tag.txt") or caExit("can't open './$asm.seqStore/readlengths-$tag.txt' for reading: $!", undef);
        while (<F>) {
            $hist .= "--    $_";
        }
        close(F);
    }

    #  Return the ASCII histogram.

    addToReport("${tag}SeqStore", $hist);
    generateReport($asm);

    #return($hist);
}



#  
#sub updateSequenceStore ($$) {
#    my $asm    = shift @_;
#    my $tag    = shift @_;
#    my $bin    = getBinDirectory();
#
#    if (! -e "$asm.seqStore/readlengths-$tag.txt") {
#        addToReport("${tag}SeqStore", generateReadLengthHistogram($tag, $asm));
#        stashSeqStore($asm);
#    }
#    generateReport($asm);
#}
