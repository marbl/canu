#!/usr/bin/perl

use strict;

my %UIDtoNAME;

#  Usage
#  fastqUIDmap fastq fastq fastq

my $fastqUIDmap = shift @ARGV;

print STDERR "Loading UID map from '$fastqUIDmap'.\n";

open (F, "< $fastqUIDmap") or die "Failed to open '$fastqUIDmap'\n";
while (<F>) {
    chomp;

    my @v = split '\s+', $_;

    if      (scalar(@v) == 3) {
        $UIDtoNAME{$v[0]} = $v[2];

    } elsif (scalar(@v) == 6) {
        $UIDtoNAME{$v[0]} = $v[2];
        $UIDtoNAME{$v[3]} = $v[5];

    } else {
        die "unknown format '$_'\n";
    }
}
close(F);

my $inFile;
my $otFile;

while (scalar(@ARGV)) {
    $inFile = shift @ARGV;
    $otFile = $inFile;

    if ($otFile =~ m/(.*).fastq/) {
        $otFile = "$1.nameFix.fastq";
    }

    open(F, "< $inFile") or die "Failed to open '$inFile' for reading\n";
    open(O, "> $otFile") or die "Failed to open '$otFile' for writing\n";

    print STDERR "Renaming '$inFile' to '$otFile'.\n";

    while (!eof(F)) {
        my $a = <F>;  chomp $a;
        my $b = <F>;
        my $c = <F>;
        my $d = <F>;

        if ($a =~ m/\@(\w+)\s/) {
            die "Didn't find UID '$1'\n"  if (!exists($UIDtoNAME{$1}));
            $a = "\@$UIDtoNAME{$1}\n";
        } else {
            print "Nope '$a'\n";
        }

        print O "$a$b$c$d";
    }

    close(F);
    close(O);
}
