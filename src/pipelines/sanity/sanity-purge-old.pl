#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

#  Purge all but the last reference assembly, and all but the last two
#  assemblies, and all but the last two that failed to complete.
#
#  Up to five assemblies are saved:
#    The reference
#    Two crashes
#    Two successful finished assemblies
#
#  No crashes are save that are earlier than the latest successful assembly.
#
#  OLDEST
#    crash
#    reference (saved)
#    crash
#    finished
#    finished
#    finished  (saved)
#    crash
#    crash
#    finished  (saved)
#    crash     (saved)
#  NEWEST

my @assemblies;
my @nightlies;

my $doPurge = ($ARGV[0] eq "purge");

#  This is automagically updated by sanity-asm-done.pl, whenever an assembly finishes.
open(F, "ls POINTERS/*last |");
while (<F>) {
    if ($_ =~ m/^POINTERS\/(.*).last$/) {
        push @assemblies, $1;
    }
}

open(F, "ls -d ????-??-??-???? |");
while (<F>) {
    chomp;
    push @nightlies, $_;
}
close(F);


foreach my $asm (@assemblies) {
    my $reference;
    my $last1;
    my $last2;
    my $crash1;
    my $crash2;

    #  Read reference pointers, find the latest.
    {
        my %ref;
        my @ref;

        $ref{"0000-00-00-0000"}++;

        if (! -e "POINTERS/$asm.reference") {
            open(F, "> POINTERS/$asm.reference");
            close(F);
        }

        open(F, "< POINTERS/$asm.reference") or die;
        while (<F>) {
            chomp;
            if (-d "$_/$asm") {
                $ref{$_}++;
            }
        }
        close(F);

        @ref = sort keys %ref;
        $reference = pop @ref;
    }

    #  Read success pointers, find the latest two.
    {
        my %ref;
        my @ref;

        $ref{"0000-00-00-0000"}++;
        $ref{"0000-00-00-0001"}++;

        if (! -e "POINTERS/$asm.last") {
            open(F, "> POINTERS/$asm.last");
            close(F);
        }

        open(F, "< POINTERS/$asm.last") or die;
        while (<F>) {
            chomp;
            if (-e "$_/$asm/$asm.qc") {
                $ref{$_}++;
            }
        }
        close(F);

        @ref = sort keys %ref;
        $last1 = pop @ref;
        $last2 = pop @ref;
    }

    #  Find the last two failures.
    {
        my %ref;
        my @ref;

        $ref{"0000-00-00-0000"}++;
        $ref{"0000-00-00-0001"}++;

        foreach my $n (@nightlies) {
            if (! -e "$n/$asm/$asm.qc") {
                $ref{$n}++;
            }
        }

        @ref = sort keys %ref;
        $crash1 = pop @ref;
        $crash2 = pop @ref;

        $crash1 = "0000-00-00-0000" if ($crash1 lt $last1);
        $crash2 = "0000-00-00-0000" if ($crash2 lt $last1);
    }

    print STDERR "REF\t$reference\tLAST\t$last1\t$last2\tCRASH\t$crash1\t$crash2\tASM\t$asm\n";

    #  Save the last TWO finished assemblies, or just one?  Comment to save two.
    $last2 = "0000-00-00-0000";

    foreach my $n (@nightlies) {
        if (-d "$n/$asm") {
            my $finished = (-e "$n/$asm/$asm.qc") ? "FINISHED" : "CRASHED";

            $finished = ($n eq $reference) ? "REFERENCE" : $finished;

            if (($n ne $reference) &&
                ($n ne $last2) &&
                ($n ne $last1) &&
                ($n ne $crash2) &&
                ($n ne $crash1)) {
                print STDERR "REMOVE\t$n/$asm\n";

                if ($doPurge) {
                    if (! -d "DEL")    { system("mkdir DEL");    }
                    if (! -d "DEL/$n") { system("mkdir DEL/$n"); }

                    #  Save some juicy bits
                    rename "$n/$asm/9-terminator/$asm.qc",   "$n/$asm.qc";
                    rename "$n/$asm/4-unitigger/$asm.cga.0", "$n/$asm.cga.0";

                    rename "$n/$asm", "DEL/$n/$asm";
                }
            } else {
                print STDERR "SAVE\t$n/$asm\n";
            }
        }
    }
}


foreach my $dir (@nightlies) {
    my $asmExist = 0;

    next if (! -d "$dir/wgs");

    foreach my $asm (@assemblies) {
        $asmExist++ if (-d "$dir/$asm");
    }

    if ($asmExist == 0) {
        print STDERR "$dir has $asmExist saved assemblies; purge source code.\n";
        if ($doPurge) {
            if (! -d "DEL")      { system("mkdir DEL");      }
            if (! -d "DEL/$dir") { system("mkdir DEL/$dir"); }

            rename "$dir/wgs", "DEL/$dir/wgs";
        }
    } else {
        print STDERR "$dir has $asmExist saved assemblies.\n";
    }
}
