#!/usr/local/bin/perl

$| = 1;

use strict;

use FindBin;
use lib "/home/walenzbp/projects/scripts";
use libBri;

my $tot = 0;
my $sma = 0;
my $big = 0;

my $smafirst = 0;
my $smalast = 0;
my $bigfirst = 0;
my $biglast = 0;

my $smaoneintronF = 0;
my $smaoneintronL = 0;
my $smaoneintronB = 0;
my $bigoneintron  = 0;

my $interiorintron = 0;

my $ff=0;
my $fc=0;
my $lf=0;
my $lc=0;

my @bigA;
my @smaA;

open(SMA,  "> sma-exon-after-big-intron");
open(BIG,  "> big-exon-after-big-intron");
open(SMAO, "> sma-exon-after-big-oneintron");
open(BIGO, "> big-exon-after-big-oneintron");

while (!eof(STDIN)) {
    $tot++;

    my %p = &libBri::readPolish(*STDIN);
    my $exonsLen = scalar(@{$p{'exons'}});
    my $firstintron = 1;

    if ($exonsLen > 1) {
        my @exons  = @{$p{'exons'}};

        my $lastC = shift @exons;
loop:
        my $thisC = shift @exons;
        my $gap   = $thisC->{'GENOMICstart'} - $lastC->{'GENOMICend'};

        if ($gap > 499999) {

            if (($firstintron) && (scalar(@exons) == 0)) {
                #  Exactly one intron
                #
                if ((($lastC->{'cDNAend'} - $lastC->{'cDNAstart'}) < 50) &&
                    (($thisC->{'cDNAend'} - $thisC->{'cDNAstart'}) < 50)) {
                    $sma++;
                    $smaoneintronB++;
                    print SMAO $p{'raw'};
                } elsif (($lastC->{'cDNAend'} - $lastC->{'cDNAstart'}) < 50) {
                    $sma++;
                    $smaoneintronF++;
                    print SMAO $p{'raw'};
                } elsif (($thisC->{'cDNAend'} - $thisC->{'cDNAstart'}) < 50) {
                    $sma++;
                    $smaoneintronL++;
                    print SMAO $p{'raw'};
                } else {
                    $big++;
                    $bigoneintron++;
                    print BIGO $p{'raw'};
                }

            } elsif ($firstintron) {
                #  First intron
                #
                if (($lastC->{'cDNAend'} - $lastC->{'cDNAstart'}) < 50) {
                    $sma++;
                    $smafirst++;
                    print SMA $p{'raw'};
                    if ($p{'matchOrientation'} eq "forward") {
                        $ff++;
                    } else {
                        $fc++;
                    }
                } else {
                    $big++;
                    $bigfirst++;
                    print BIG $p{'raw'};
                }
            } elsif (scalar(@exons) == 0) {
                #  Last intron
                #
                if (($thisC->{'cDNAend'} - $thisC->{'cDNAstart'}) < 50) {
                    $sma++;
                    $smalast++;
                    print SMA $p{'raw'};
                    if ($p{'matchOrientation'} eq "forward") {
                        $lf++;
                    } else {
                        $lc++;
                    }
                } else {
                    $big++;
                    $biglast++;
                    print BIG $p{'raw'};
                }
            } else {
                #  Interior intron
                #
                $interiorintron++;
            }

            print "int: $interiorintron sma: $sma(First:$smafirst,Last:$smalast,";
            print "oneF:$smaoneintronF,oneL:$smaoneintronL,oneB:$smaoneintronB) -- big ";
            print "$big(First:$bigfirst,Last:$biglast,One:$bigoneintron) -- tot $tot -- ";
            print "ff=$ff,fc=$fc lf=$lf,lc=$lc\n";
        }

        $firstintron = 0;

        $lastC = $thisC;

        goto loop if (scalar(@exons) > 0);
    }

}

close(SMA);
close(BIG);
close(SMAO);
close(BIGO);
