#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

package sim4polish;

use strict;
use POSIX "sys_wait_h";

$| = 1;

sub import () {
}


######################################################################
#
#  Returns a modified 'raw' string, using the current values for the
#  info line.  DOES NOT rewrite the exons.
#
sub updatePolishInfoLine {
    my %p = @_;
    my @L = split '\n', $p{'raw'};
    my $l;

    shift @L;
    shift @L;

    $l  = "sim4begin\n";
    $l .= "$p{'estID'}\[$p{'estLen'}-$p{'pA'}-$p{'pT'}\] ";
    $l .= "$p{'dbID'}\[$p{'dbLo'}-$p{'dbHi'}\] ";
    $l .= "<$p{'numMatches'}-$p{'numMatchesN'}-$p{'percentID'}-$p{'matchOrientation'}-$p{'strandPrediction'}>\n";

    foreach my $x (@L) {
        $l .= "$x\n";
    }

    return($l);
}


######################################################################
#
#  Subroutine to read a single sim4 polish, and return it as a structure.
#
sub readPolish (*) {
    local *READPOLISHFH = shift;
    my %p;
    my $line;
    my $save;

    #  These are the fields returned.
    #
    $p{'raw'} = undef;

    $p{'estID'} = undef;
    $p{'estDefLine'} = undef;
    $p{'estLen'} = undef;
    $p{'pA'} = undef;
    $p{'pT'} = undef;

    $p{'dbID'} = undef;
    $p{'dbDefLine'} = undef;
    $p{'dbLen'} = undef;
    $p{'dbLo'} = undef;
    $p{'dbHi'} = undef;

    $p{'comment'} = undef;

    $p{'numMatches'} = undef;
    $p{'numMatchesN'} = undef;
    $p{'percentID'} = undef;
    $p{'coveragebp'} = undef;
    $p{'coverage'} = undef;
    $p{'matchOrientation'} = undef;
    $p{'strandPrediction'} = undef;

    #  An array of references to hashes, one hash for each exon.
    $p{'exons'} = ();


    #  Skip lines until the next match.  If used properly, on a proper
    #  file, this should be silent.  After the loop, we are at the
    #  start of a polish; the line should be "sim4begin".
    #
    $line = <READPOLISHFH>;
    while (defined($line) && ($line !~ m/^sim4begin$/)) {
        chomp $line;
        print STDERR "Skipped: '$line'\n";
        $line = <READPOLISHFH>;
    }
    $save = $line;

    #  Return now if were are out of file
    #
    return(%p) if (eof(READPOLISHFH));


    #  Read the description line
    #
    $line  = <READPOLISHFH>;
    $save .= $line;

    if ($line =~ m/^(\d+)\[(\d+)-+(\d+)-+(\d+)\]\s+(\d+)\[(\d+)-(\d+)\]\s+\<(\d+)-(\d+)-(\d+)-(\w+)-(\w+)\>$/) {
        $p{'estID'}            = $1;
        $p{'estLen'}           = $2;
        $p{'pA'}               = $3;
        $p{'pT'}               = $4;
        $p{'dbID'}             = $5;
        $p{'dbLo'}             = $6;
        $p{'dbHi'}             = $7;
        $p{'numMatches'}       = $8;
        $p{'numMatchesN'}      = $9;
        $p{'percentID'}        = $10;
        $p{'matchOrientation'} = $11;
        $p{'strandPrediction'} = $12;
    } else {
        print STDERR "expecting description line, got: '$line'\n";
        return(%p);
    }


    #  Read the two deflines, if they exist.
    #
    $line = <READPOLISHFH>;

    if ($line =~ m/^comment=\s*(.*)\s*$/) {
        $p{'comment'} = $1;
        $save .= $line;
        $line  = <READPOLISHFH>;
    } else {
        #print STDERR "libBri::readPolish()-- WARNING:  Didn't get comment!\n";
        #print STDERR "libBri::readPolish()-- WARNING:  $line";
    }
    if ($line =~ m/^edef=(.*)$/) {
        $p{'estDefLine'} = $1;
        $save .= $line;
        $line  = <READPOLISHFH>;
    } else {
        #print STDERR "libBri::readPolish()-- WARNING:  Didn't get edef!\n";
        #print STDERR "libBri::readPolish()-- WARNING:  $line";
    }

    if ($line =~ m/^ddef=(.*)$/) {
        $p{'dbDefLine'} = $1;
        $save .= $line;
        $line  = <READPOLISHFH>;
    } else {
        #print STDERR "libBri::readPolish()-- WARNING:  Didn't get ddef!\n";
        #print STDERR "libBri::readPolish()-- WARNING:  $line";
    }


    #  Read the exons
    #
    my $exonAlign      = 0;
    my $exonAlignFirst = 1;
    my $exonCoverage   = 0;

    while (defined($line) && ($line !~ m/^sim4end$/)) {

        #  If this match succeeds, we have an exon description.
        #  Otherwise, it's an alignment line.
        #
        if ($line =~ /^(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+\<(\d+)-(\d+)-(\d+)\>\s+(.*)$/) {
            my $e = {};

            $exonCoverage += $2 - $1 + 1;

            $e->{'cDNAstart'}         = $1;
            $e->{'cDNAend'}           = $2;
            $e->{'GENOMICstart'}      = $3;
            $e->{'GENOMICend'}        = $4;
            $e->{'numMatches'}        = $5;
            $e->{'numMatchesN'}       = $6;
            $e->{'percentID'}         = $7;
            $e->{'intronOrientation'} = $8;

            push @{$p{'exons'}}, $e;
        } else {
            if ($exonAlignFirst) {
                $p{'exons'}[$exonAlign]->{'cDNAalign'} = $line;
                chomp $p{'exons'}[$exonAlign]->{'cDNAalign'};
                $exonAlignFirst = 0;
            } else {
                $p{'exons'}[$exonAlign]->{'GENOMICalign'} = $line;
                chomp $p{'exons'}[$exonAlign]->{'GENOMICalign'};
                $exonAlignFirst = 1;
                $exonAlign++;
            }
        }

        $save .= $line;
        $line  = <READPOLISHFH>;
    }

    $save .= $line;

    if (($p{'pA'} + $p{'pT'}) >= $p{'estLen'}) {
        $p{'coverage'} = 0;
    } else {
        $p{'coveragebp'} = $exonCoverage;
        $p{'coverage'}   = 100.0 * $exonCoverage / ($p{'estLen'} - $p{'pA'} - $p{'pT'});
    }

    $p{'raw'} = $save;

    return(%p);
}



1;
