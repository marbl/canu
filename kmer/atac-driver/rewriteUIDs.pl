#!/usr/bin/perl

#  Reads an atac file with atac-format IDs, writes an atac file with
#  UIDs (the first word in the defline).  This is the last step in the
#  normal atac pipeline.

use strict;

my $atacfile = shift @ARGV;

my $seqA;
my $tagA;
my %uidA;

my $seqB;
my $tagB;
my %uidB;

my $iid;

open(F, "< $atacfile") or die;
while (!defined($seqA) || !defined($tagA) || !defined($seqB) || !defined($tagB)) {
    $_ = <F>;
    $seqA = $1 if (m/^\/assemblyFile1=(.*)$/);
    $tagA = $1 if (m/^\/assemblyId1=(.*)$/);
    $seqB = $1 if (m/^\/assemblyFile2=(.*)$/);
    $tagB = $1 if (m/^\/assemblyId2=(.*)$/);
}
close(F);

if (!defined($seqA) || !defined($tagA) || !defined($seqB) || !defined($tagB)) {
    die "Something fishy.  Didn't find seqs or tags in '$atacfile'.\n";
}

$iid = 0;
open(F, "< $seqA") or die "Failed to open '$seqA'\n";
while (<F>) {
    if (m/^>(\S+)\s*.*$/) {
        #chomp;
        #print STDERR "$tagA:$iid -> $_\n";
        $uidA{"$tagA:$iid"} = $1;
        $iid++;
    }
}
close(F);


$iid = 0;
open(F, "< $seqB") or die "Failed to open '$seqA'\n";
while (<F>) {
    if (m/^>(\S+)\s*.*$/) {
        #chomp;
        #print STDERR "$tagB:$iid -> $_\n";
        $uidB{"$tagB:$iid"} = $1;
        $iid++;
    }
}
close(F);


$, = " ";
$\ = "\n";

open(F, "< $atacfile") or die;
while (<F>) {
    chomp $_;

    my @v = split '\s+', $_;

    if (m/^M/) {
        die "Didn't find uidA for $v[4]\n" if (!defined($uidA{$v[4]}));
        die "Didn't find uidB for $v[8]\n" if (!defined($uidB{$v[8]}));

        $v[4] = $uidA{$v[4]};
        $v[8] = $uidB{$v[8]};
        print @v;
    } else {
        print $_;
    }
}
close(F);
