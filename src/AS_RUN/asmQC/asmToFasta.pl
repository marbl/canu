#!/usr/local/bin/perl

use strict;
use XML::Parser;
use Getopt::Long;
use TIGR::AsmLib;

my ($orgID,$orgXML,$useScafOrient,$gcode,$pref);
my $degenOut = '';
my $usage = $0 .' [ARGS] < assembly.asm > assembly.fasta
ARGS are all required:
--xml Moore130.xml
--indexInXml M006
--useScafOrient [true|false]
--prefix output_file_prefix [M006 will make M006.contigs.qv M006.contigs.fasta]
[--gcode 1 for most eucaryotes, default is 11 for bacteria ]
[--degen output degenerates instead of scaffold contigs ]

If --useScafOrient is false the orientation in the scaffold will be ignored, and
the contigs will be output as they appear in the assembly file.
If --useScafOrient is true the contigs will be reverse complemented before being
output if the scaffold structure shows them in the reverse orientation.

Sample XML File:
<?xml version="1.0" encoding="ISO-8859-1"?>

<table>
<column name="Library Prefix">Macaca</column>
<column name="Strain"></column>
<column name="Species">Macaca Mulatta</column>
</table>
';

die $usage unless @ARGV;

GetOptions('xml|x=s'=> \$orgXML,'indexInXml|i=s'=> \$orgID,'prefix|p=s'=>\$pref,
    'useScafOrient|u=s' => \$useScafOrient, 'gcode|g=i' => \$gcode,
    'degen|d' => \$degenOut
) || die $usage ;

if ($useScafOrient eq 'true' ) {
    $useScafOrient = 1;
} elsif ($useScafOrient eq 'false') {
    $useScafOrient = 0;
} else {
    die $usage;
}
$gcode = 11 unless defined $gcode;
die $usage unless defined $pref;
open(CONTIGS,">$pref.contigs.fasta") || die "Can't write to $pref.contigs.fasta";
open(QUALS,">$pref.contigs.qv") || die "Can't write to $pref.contigs.qv";

my %contigs;
my %quals;
my %seen;
my $defLine;
my $offset=0;
my $qvPerLine = 25;
my $minGapSize = 20;
sub printContig($$;$) {
    my ($id,$len,$scafID) = @_;
    my $scaf = '';
    $scaf = " [scaffold=$scafID]" if defined $scafID;
    my $def = $defLine."$scaf [offset=$offset] [len=$len]";
    if ( ! $seen{$id} ) {
        printFastaSequence(\*CONTIGS,$id.$def, $contigs{$id});
        my $off = 0;
        my $baseLeft = $len > $qvPerLine ? $qvPerLine : $len ;
        print QUALS ">$id $def\n";
        while( my $qv60 = substr $quals{$id}, $off, $baseLeft) {
            for(my $i=0;$i < $baseLeft; $i++) {
                my $intQv = ord(substr $qv60,$i,1) - ord("0");
                die "$i $intQv $baseLeft $qv60." if $intQv < 0;
                printf QUALS "%02d ",$intQv;
            }
            print QUALS "\n";
            $off += $baseLeft;
            $baseLeft = $len - $off > $qvPerLine ? $qvPerLine : $len - $off ;
        }
        $seen{$id} = 1;
    }
}
sub outputContigs($$$$;$) {
    my ($id1,$id2,$mean,$orient,$scafID) = @_;
    my $gap = $mean > $minGapSize ? int($mean) : $minGapSize;

    if ($orient eq "A") {
        $contigs{$id1} = reverseComplement($contigs{$id1});
        $contigs{$id2} = reverseComplement($contigs{$id2});
        $quals{$id1}   = reverse $quals{$id1};
        $quals{$id2}   = reverse $quals{$id2};
    } elsif ($orient eq "I") {
        $contigs{$id2} = reverseComplement($contigs{$id2});
        $quals{$id2}   = reverse $quals{$id2};
    } elsif ($orient eq "N") {
    } elsif ($orient eq "O") {
        $contigs{$id1} = reverseComplement($contigs{$id1});
        $quals{$id1}   = reverse $quals{$id1};
    } else {
        die "Unknown orientation $orient";
    }
    my $len = length($contigs{$id1});
    printContig($id1,$len,$scafID);

    $offset += $len + $gap;

    $len = length($contigs{$id2});
    printContig($id2,$len,$scafID);
}
my %cols;
my $col;
my $found=0;
sub cols{
    my $p = shift;
    $cols{$col} = shift;
    $found = 1 if $col eq 'Library Prefix' && $cols{$col} eq $orgID;
    $p->setHandlers(Char => undef);
}
sub startXML{
   my $p = shift;
   if ($_[0] eq 'column') {
       my $c = $_[2];
       if ($c eq 'Strain' || $c eq 'Species' || $c eq 'Library Prefix') {
           $col = $c;
           $p->setHandlers(Char => \&cols);
       } elsif ($c eq '16s rRNA' && $found) {
           my ($species,$strain) = ($cols{Species},$cols{Strain});
           my $match = "\$species =~ s'$strain''";
           eval "$match";
           $species =~ s/\s+$//;
           $defLine = " [organism=$species] [strain=$strain] [gcode=$gcode]";

           $p->setHandlers(Start => undef);
       }
   }
}
# Begin the main program
my $xmlP = new XML::Parser(Handlers => {Start=>\&startXML});
$xmlP->parsefile($orgXML);

while (my $record = getCARecord(\*STDIN)){
    my ($rec, $fields, $recs) = parseCARecord($record);
    if ($rec eq "CCO"){
        my $id = getCAId($$fields{acc});

        my $len = $$fields{len};
        my $seq = $$fields{cns};
        my $qlt = $$fields{qlt};
        my $idx = index $seq,'-';
        my @dashPos;
        while ( $idx > 0 ) {
            push @dashPos,$idx;
            $idx = index $seq,'-',$idx+1;
        }
        for $idx (reverse @dashPos) {
            substr $qlt,$idx,1,'';
        }
        $seq =~ s/[-\s]//g;
        $qlt =~ s/\s//g;
        $contigs{$id} = $seq;
        $quals{$id} = $qlt;

    } elsif (!$degenOut && $rec eq "SCF") {
        $offset = 0;
        $$fields{acc} =~ /\((\S+),\S+/;
        $$fields{acc} =~ /(\S+)/ unless defined $1 && $1 ne '';
        my $scfId = $1;
        for my $sr (@$recs) {
            my ($sub,$sfs,$srecs) = parseCARecord($sr);
            #print STDERR "$sub $sfs->{ct1} $sfs->{ct2} $sfs->{mea} $sfs->{std} $sfs->{ori}\n";
            my $ori = 'N';
            $ori = $sfs->{ori} if $useScafOrient;
            outputContigs($sfs->{ct1},$sfs->{ct2},$sfs->{mea},$ori,$scfId);
        }
    } elsif ($degenOut && $rec eq 'DSC') {
        $offset = 0;
        $$fields{acc} =~ /\((\S+),\S+/;
        my $scfId = $1;
        outputContigs($fields->{ctg},$fields->{ctg},0,'N',$scfId);
    }

}
# don't include unplaced contigs
#for my $id (keys %contigs) {
#    print(">$id\n",$contigs{$id},"\n") unless $seen{$id};
#}
