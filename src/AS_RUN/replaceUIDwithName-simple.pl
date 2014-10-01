#!/usr/bin/perl

use strict;

my %UIDtoNAME;

#  Usage
#  fastqUIDmap posmnap posmap posmap

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

    print STDERR "Renaming '$inFile' to '$inFile.UID'.\n";

    rename "$inFile", "$inFile.UID";

    open(F, "< $inFile.UID") or die "Failed to open '$inFile.UID' for reading\n";
    open(O, "> $inFile")     or die "Failed to open '$inFile' for writing\n";

    while (!eof(F)) {
        my $a = <F>;
        my @a = split '\s+', $a;

        foreach my $a (@a) {
            if (exists($UIDtoNAME{$a})) {
                $a = $UIDtoNAME{$a};
            }
        }

        print O join "\t", @a;
        print O "\n";
    }

    close(F);
    close(O);
}
