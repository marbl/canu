#!/usr/local/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use TIGR::AsmLib;

=head 1

This program produces an agp file needed for WGS submission to genbank.
It does not yet add a gap at the end for a circular genome.

=cut

my $usage = "$0 [-d] < Assembly.asm > assembly.agp
The -d option causes degenerates to be output instead of normal scaffolds.";

die $usage if @ARGV > 1 || @ARGV == 1 && $ARGV[0] ne '-d';

my $outputDegen = 0;
$outputDegen = 1 if @ARGV == 1;

my %contigs;
my %seen;
my $defLine;
my $offset=0;
my $cnt=0;

sub printContig($$$$) {
    my ($scf,$id,$len,$ort) = @_;
    if ( ! $seen{$id} ) {
        print $scf,$offset+1,$offset+$len,++$cnt,'W',$id,1,$len,$ort;
        $seen{$id} = 1;
    }
}
sub outputContigs($$$$$) {
    my ($scf,$id1,$id2,$mean,$orient) = @_;
    # set unknown and small gaps to size 100, since that's what NCBI asked for
    my $gap = $mean > 20 ? int($mean) : 20;

    my ($ort1,$ort2) = ('+','+');
    if ($orient eq "A") {
        ($ort1,$ort2) = ('-','-');
    } elsif ($orient eq "I") {
        $ort2 = '-';
    } elsif ($orient eq "N") {
    } elsif ($orient eq "O") {
        $ort1 = '-';
    } else {
        die "Unknown orientation $orient";
    }
    local $, = "\t";
    local $\ = "\n";
    my $len = $contigs{$id1};
    printContig($scf,$id1,$len,$ort1);

    $offset += $len;

    print $scf,$offset+1,$offset+$gap,++$cnt,'N',$gap,'fragment','yes' if $id1 ne $id2;

    $offset += $gap;

    $len = $contigs{$id2};
    printContig($scf,$id2,$len,$ort2);
}

while (my $record = getCARecord(\*STDIN)){
    my ($rec, $fields, $recs) = parseCARecord($record);
    if ($rec eq 'CCO'){
        my $id = getCAId($$fields{acc});

        my $len = $$fields{len};
        my $seq = $$fields{cns};
        $seq =~ s/[-\s]//g;
        $contigs{$id} = length($seq);

    } elsif (!$outputDegen && $rec eq 'SCF') {
        my $scfID = getCAId($$fields{acc});
        $offset = 0;
        $cnt = 0;
        for my $sr (@$recs) {
            my ($sub,$sfs,$srecs) = parseCARecord($sr);
            #print STDERR "$sub $sfs->{ct1} $sfs->{ct2} $sfs->{mea} $sfs->{std} $sfs->{ori}\n";
            outputContigs($scfID,$sfs->{ct1},$sfs->{ct2},$sfs->{mea},$sfs->{ori});
        }
    } elsif ($outputDegen && $rec eq 'DSC') {
        my $scfID = getCAId($$fields{acc});
        $offset = 0;
        $cnt = 0;
        outputContigs($scfID,$fields->{ctg},$fields->{ctg},0,'N');
    }

}
