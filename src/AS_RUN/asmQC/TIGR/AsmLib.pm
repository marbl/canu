# $Id: AsmLib.pm,v 1.1 2006-09-04 07:15:37 brianwalenz Exp $
#
# File: AsmLib.pm
# Authors: Mihai Pop
#
#  Copyright @ 2002, The Institute for Genomic Research (TIGR).  All
#  rights reserved.
#
#  This software is provided "AS IS".  TIGR makes no warranties, express
#  or implied, including no representation or warranty with respect to
#  the performance of the software and derivatives or their safety,
#  effectiveness, or commercial viability.  TIGR does not warrant the
#  merchantability or fitness of the software and derivatives for any
#  particular purpose, or that they may be exploited without infringing
#  the copyrights, patent rights or property rights of others. TIGR shall
#  not be liable for any claim, demand or action for any loss, harm,
#  illness or other damage or injury arising from access to or use of the
#  software or associated information, including without limitation any
#  direct, indirect, incidental, exemplary, special or consequential
#  damages.

#  This software program may not be sold, leased, transferred, exported
#  or otherwise disclaimed to anyone, in whole or in part, without the
#  prior written consent of TIGR.
#
#
# Routines that help processing Fasta and Celera .asm files.
#

=head1 NAME

AsmLib - A library of perl routines for processing Celera .frg, .asm files and 
TIGR .contig, .seq, and .qual files.

=head1 SYNOPSIS

use TIGR::AsmLib;

=head1 DESCRIPTION

A set of functions that help the processing of Celera .asm files and
Fasta files as well as the creation of TIGR .asm and .contig files.

=cut

package TIGR::AsmLib;

use strict;

## configuration management            
our $VERSION = "2.00";
our $VERSION_STRING = "$VERSION (Build " . (qw/$Revision: 1.1 $/ )[1] . ")";
our @DEPEND = ();

BEGIN {
    use Exporter ();
    use vars qw(@EXPORT @EXPORT_OK @ISA %EXPORT_TAGS);

    @ISA         = qw(Exporter);
    @EXPORT      = qw(&getCARecord
                      &parseCARecord
                      &getCAId
                      &printContigRecord
                      &printSequenceRecord
                      &printFastaSequence
                      &printFastaQual
                      &reverseComplement
                      $VERSION
                      $VERSION_STRING
                      @DEPEND);
    %EXPORT_TAGS = ();
    @EXPORT_OK   = ();
}

use vars @EXPORT;
use vars @EXPORT_OK;


##############################################

=over

=item B<my $rec = getCARecord(\*STDIN);>

Reads from stdin the text between "extreme" { and } .

 for example:

  {A
    {B
    }
   }
 
Returns the whole: {A{B}}

=cut

sub getCARecord
{
    my $file = shift;
    
    my $level = 0;
    my $block = "";

    local $/ = "\n{";
    while (<$file>){
        $level++;
        my $closeCnt = substr($_,-7) =~ tr/}/}/;
        $level -= $closeCnt;

        $block .= $_;
#        push @block,$_;
        if ($level == 0){
            chomp $block;
            last;
        }
    }

    if ($level != 0){
      die ("end of file reached before end of block\n");
    }

    if ($block ne ""){
#      return join '',@block;
      return $block;
    } else {
      return undef;
    }
}# getCARecord

sub print_sequence
{
    my $file = shift;
    my $seqs = shift;

    for (my $j = 0; $j < length($seqs); $j += 60){
        print $file substr($seqs, $j, 60), "\n";
    }
}


######################################################3

=item B<my($id, $fields, $recs) = parseCARecord($rec);>

Parses a record and returns a triplet consisting of
   - record type
   - hash of fields and values
   - array of sub-records

=cut

sub parseCARecord
{
    my $record = shift;

    my @lines = split('\n', $record);

    my $type;
    my %fields;
    my @recs;

    # get record type
    $lines[0] =~ /^\{?(\w+)/;
    if (! defined $1){
        die ("Wierd start of record: $record\n");
    }
    $type = $1;

    if ($lines[$#lines] !~ /^\s*\}/){
      die ("Wierd end of record: $record\n");
    }

    my $level = 0;
    my $fieldname;
    for (my $i = 1; $i < $#lines; $i++){
      if ($lines[$i] =~ /^(\w+):(\S+)$/){   # simple field
          $fields{$1} = $2;
      } # simple field
      elsif ($lines[$i] =~ /^(\w+):$/){ # complex field
          $fieldname = $1;
          $fields{$fieldname} = "";
          $i++;
          while ($i < $#lines && ($lines[$i] !~ /^\.$/)){
            $fields{$fieldname} .= "$lines[$i]\n";
            $i++;
          }
      } # complex field
      elsif ($lines[$i] =~ /^\s*\{/){ # subrecord
          $level = 1;

          my $thisrec = ++$#recs;
          
          $recs[$thisrec] = "$lines[$i]\n";
          $i++;
          while ($level > 0 && $i < $#lines){
            if ($lines[$i] =~ /^\s*\{/){
                $level++;
            }
            if ($lines[$i] =~ /^\s*\}/){
                $level--;
            }
            $recs[$thisrec] .= "$lines[$i]\n";
            if ($level == 0){
                last;
            } else {
                $i++;
            }
          }

          if ($level != 0){
            die ("Error parsing sub_record in: $record\n");
          }
      } # subrecord
    } # for $i...
    
    return ($type, \%fields, \@recs);
} # parseCARecord

################################################

=item B<my($id) = getCAId($CAid);>

Obtains the ID from a "paired" id, that is, converts (10, 1000) into 10.
If the Id is not a pair in parantheses, it returns the input.
Thus, getCAId('(10, 1000)') returns 10 while getCAId("abba") returns "abba".

=cut

sub getCAId
{
    my $string = shift;

    if (defined $string && $string =~ /\((\d+),(\d+)\)/){
      return $1;
    } else {
      return $string; # just in case we have a real ID
    }
} # getCAId

################################################

=item B<printContigRecord($file, $id, $len, $nseq, $sequence, $how);>

Prints contig in specified format
Inputs are:
   $file - output file (opened for writing)
   $id  - contig ID
   $len - contig length
   $nseq - number of sequences in contig (same as number of sequence records
that will follow the contig
   $sequence - consensus sequence for the contig
   $how - what type of output is required:
        contig - TIGR .contig format
        asm    - TIGR .asm format
        fasta  - multi-fasta format

=cut

sub printContigRecord
{
    my $file = shift;
    my $id = shift;
    my $len = shift;
    my $nseq = shift;
    my $sequence = shift;
    my $how = shift;

    if ($how eq "contig"){
      print $file "\#\#$id $nseq $len bases, 00000000 checksum.\n";
      print_sequence($file, $sequence);
      return;
    } # if $how eq "contig"
    elsif ($how eq "asm"){
      my $strip = $sequence;
      $strip =~ s/-//g;
# get the current date
      my $date = `date +'%D %T'`;
      chomp $date;

      my $quality = "0x";
      for (my $i = 0; $i < length($sequence); $i++){
          $quality .= "06";
      }

      print $file "sequence\t$strip\n";
      print $file "lsequence\t$sequence\n";
      print $file "quality\t$quality\n";
      print $file "asmbl_id\t$id\n";
      print $file "seq_id\t\n";
      print $file "com_name\t\n";
      print $file "type\t\n";
      print $file "method\tCelera Assembler\n";
      print $file "ed_status\t\n";
      print $file "redundancy\t\n";
      print $file "perc_N\t\n";
      print $file "seq\#\t$nseq\n";
      print $file "full_cds\t\n";
      print $file "cds_start\t\n";
      print $file "cds_end\t\n";
      print $file "ed_pn\t$ENV{USER}\@$ENV{HOSTNAME}\n";
      print $file "ed_date\t$date\n";
      print $file "comment\t\n";
      print $file "frameshift\t\n";
      return;
    } # if $how eq "asm"
    elsif ($how eq "fasta")
    {
      print $file ">$id $nseq $len bases\n";
      print_sequence($file, $sequence);
      return;
    }
} # printContigRecord

################################################

=item B<printSequenceRecord($file, $name, $seq, $offset, $rc, $seqleft, $seqright, $asml, $asmr, $type);>

Prints the record for a sequence aligned to a contig
Inputs are:
   $file - output file opened for writing
   $name - sequence name
   $seq - actual sequence
   $offset - offset in consensus
   $rc - "RC" if sequence is reverse complemented, "" otherwise
   $seqleft, $seqright - alignment range within sequence
   $asml, $asmr - alignment range within consensus
   $type - type of output:
           contig -  output is in TIGR .contig format 
           asm    -  output is in TIGR .asm format

=cut

sub printSequenceRecord
{
    my($file, $name, $seq, $offset, $rc, 
       $seqleft, $seqright, $asml, $asmr, $type) = @_;

    if ($type eq "contig"){
      print $file "\#$name($offset) [$rc] ", 
      length($seq), 
      " bases, 00000000 checksum. {$seqleft $seqright} <$asml $asmr>\n";
      
      print_sequence($file, $seq);
    }

    if ($type eq "asm"){
      print $file "\n";
      print $file "seq_name\t$name\n";
      print $file "asm_lend\t$asml\n";
      print $file "asm_rend\t$asmr\n";
      print $file "seq_lend\t$seqleft\n";
      print $file "seq_rend\t$seqright\n";
      print $file "best\t\n";
      print $file "comment\t\n";
      print $file "db\t\n";
      print $file "offset\t$offset\n";
      print $file "lsequence\t$seq\n";
    }

    return;
} # printSequenceRecord

################################################

=item B<printFastaSequence($file, $header, $seq);>

Prints sequence in Fasta format
Inputs are:
   $file - output file opened for writing
   $header - Fasta header (without >)
   $seq - sequence to be written

=cut

sub printFastaSequence
{
    my($file) = $_[0];
    my($header) = $_[1];
    my($seqs) = $_[2];

    print $file ">$header\n";
    print_sequence($file, $seqs);
} # printFastaSequence

################################################

=item B<printFastaQual($file, $header, $qual);>

Prints quality values in Fasta format.
Inputs are:
    $file - output file
    $header - fasta header (without >)
    $qual - string of quality values

=cut

sub printFastaQual
{
    my($file) = $_[0];
    my($header) = $_[1];
    my($quals) = $_[2];
    my(@qv);
 
    print $file ">$header\n";

    @qv = split(' ', $quals);
    
    for (my $j = 0; $j <= $#qv; $j += 17){
        print $file join(" ", @qv[$j .. $j + 16]), "\n";
    }
} # printFastaQual


################################################

=item B<my($rev) = reverseComplement($seq);>

Reverse complements a sequence.

=cut

sub reverseComplement($) {
    my($string) = @_;

    $string = reverse($string);
    $string =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBbXxNn.-/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVvXxNn.-/;
    return $string;
} # reverseComplement

=back

=cut

END{}

1;
