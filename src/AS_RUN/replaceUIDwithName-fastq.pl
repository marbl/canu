#!/usr/bin/perl

use strict;

#  Usage
#  fastqUIDmap fastq fastq fastq

my $fastqUIDmap = shift @ARGV;

my $inFile;
my $otFile;

my %UIDtoNAME;

my $namesLoaded   = 0;
my $readsRenamed  = 0;

sub loadMoreNames {
    my $lim = 100000;

    undef %UIDtoNAME;

    while (<N>) {
        chomp;

        my @v = split '\s+', $_;

        if      (scalar(@v) == 3) {
            $UIDtoNAME{$v[0]} = $v[2];
            $namesLoaded++;

        } elsif (scalar(@v) == 6) {
            $UIDtoNAME{$v[0]} = $v[2];
            $UIDtoNAME{$v[3]} = $v[5];
            $namesLoaded++;
            $namesLoaded++;

        } else {
            die "unknown format '$_'\n";
        }

        if (--$lim == 0) {
            return;
        }
    }
}



while (scalar(@ARGV)) {
    $inFile = shift @ARGV;
    $otFile = $inFile;

    if ($inFile =~ m/(.*).fastq/) {
        $otFile = "$1.fastq.RENAMING";
    } else {
        die "Failed to generate output file name.\n";
    }

    open(F, "< $inFile")      or die "Failed to open '$inFile' for reading\n";
    open(O, "> $otFile")      or die "Failed to open '$otFile' for writing\n";
    open(N, "< $fastqUIDmap") or die "Failed to open '$fastqUIDmap'\n";

    $namesLoaded   = 0;
    $readsRenamed  = 0;

    print STDERR "Renaming '$inFile' to '$otFile'.\n";

    while (!eof(F)) {
        my $a = <F>;  chomp $a;
        my $b = <F>;
        my $c = <F>;
        my $d = <F>;

        if ($a =~ m/\@(\w+),\w+\s*/) {
            #  UID,IID
            $a = $1;

        } elsif ($a =~ m/\@(\w+)\s*/) {
            #  UID
            $a = $1;

        } else {
            die "Nope '$a'\n";
        }

        while (!exists($UIDtoNAME{$a})) {
            loadMoreNames();

            #if ((!exists($UIDtoNAME{$a})) && ($readsRenamed > 0)) {
            #    print STDERR "WARNING:  Looping to load more names; out of sync?\n";
            #}
        }

        die "Didn't find UID '$a'\n"  if (!exists($UIDtoNAME{$a}));
        $a = "\@$UIDtoNAME{$a}\n";
        
        print O "$a$b$c$d";

        $readsRenamed++;

        if (($readsRenamed % 10000) == 0) {
            print STDERR "Renamed $readsRenamed reads using $namesLoaded names.\r";
        }
    }

    print STDERR "Renamed $readsRenamed reads using $namesLoaded names.\n";

    close(F);
    close(O);
    close(N);

    if ($inFile =~ m/(.*).fastq/) {
        rename "$1.fastq",          "$1.CA_UIDs.fastq";
        rename "$1.fastq.RENAMING", "$1.fastq";
    }
}
