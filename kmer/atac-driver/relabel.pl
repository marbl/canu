#!/usr/bin/perl

#  Reads an atac file, relabels the sequence names (e.g., WGSA:4) with
#  the defline ID's.

use strict;

sub readDeflines ($) {
    my $file = $_[0];
    if ($file =~ m/^(.*).fasta/) {
        $file = $1;
    }

    if (-e "$file.deflines") {
        $file = "$file.deflines";
    } elsif (-e "$file.fasta.deflines") {
        $file = "$file.fasta.deflines";
    } else {
        print STDERR "Dang, gotta grep the deflines!\n";
        system("grep '>' $file.fasta > $file.deflines");
        $file = "$file.deflines";
    }

    my @nameA;

    #print STDERR "$file\n";

    open(Z, "< $file") or die "Failed to open '$file'\n";
    while (!eof(Z)) {
        my $n = <Z>;
        if ($n =~ m/^\>\s*(\S+)\s*/) {
            push @nameA, $1;
        } else {
            chomp $n;
            print STDERR "Failed to match defline '$n'\n";
        }
    }
    close(Z);

    return(@nameA);
}


my $file = shift @ARGV;
my @nameA;
my @nameB;

if ($file eq "-A") {
    @nameA = readDeflines(shift @ARGV);
    $file = shift @ARGV;
}
if ($file eq "-B") {
    @nameB = readDeflines(shift @ARGV);
    $file = shift @ARGV;
}

open(F, "< $file")      or die "Failed to open '$file' for input\n";
open(G, "> $file.uids") or die "Failed to open '$file.uids' for output\n";

while (<F>) {
    if (m/assemblyFile1=(.*)$/) {
        @nameA = readDeflines($1);
    print STDERR "num nameA = ", scalar(@nameA), "\n";
    }
    if (m/assemblyFile2=(.*)$/) {
        @nameB = readDeflines($1);
    print STDERR "num nameB = ", scalar(@nameB), "\n";
    }

    if (m/^M\s/) {
        my @v = split '\s+', $_;

        if ($v[4] =~ m/^\w+:(\d+)$/) {
            if (defined($nameA[$1])) {
                $v[4] = $nameA[$1];
            } else {
                die "Didn't find nameA for $1\n";
            }
        } else {
            die "Didn't match v[4] = $v[4]\n";
        }

        if ($v[8] =~ m/^\w+:(\d+)$/) {
            if (defined($nameB[$1])) {
                $v[8] = $nameB[$1];
            } else {
                die "Didn't find nameA for $1\n";
            }
        } else {
            die "Didn't match v[8] = $v[8]\n";
        }

        #  Special case stuff....
        #
        if      ($v[4] =~ m/^Chr(\d+)$/) {
            $v[4] = "mchr$1";
        } elsif ($v[4] =~ m/^Chr(\d+)_random$/) {
            $v[4] = "mchr${1}r";
        } elsif ($v[4] =~ m/^SCAFFOLD(\d+)$/) {
            $v[4] = "bscf$1";
        } elsif ($v[4] =~ m/^Contig(\d+)$/) {
            $v[4] = "wscf$1";
        } elsif ($v[4] =~ m/^chr(\d+)$/) {
            $v[4] = "hchr$1";
        }

        if      ($v[8] =~ m/^SCAFFOLD(\d+)$/) {
            $v[8] = "bscf$1";
        } elsif ($v[8] =~ m/^Contig(\d+)$/) {
            $v[8] = "wscf$1";
        } elsif ($v[8] =~ m/^chr(\d+)$/) {
            $v[8] = "hchr$1";
        }



        my $line = join " ", @v;
        print G "$line\n";
    } else {
        print G $_;
    }
}

close(G);
close(F);
