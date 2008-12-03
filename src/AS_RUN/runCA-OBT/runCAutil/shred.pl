#!/usr/bin/perl

use strict;
use FileHandle;

#
#  Parameters
#

my $MIN_COVERAGE      = 1;  #  Should be 2 if there are "fake" reads in ace file

my $MIN_READS         = 4;
my $MIN_CONTIG_SIZE   = 600;

my $SHRED_READ_LENGTH = 600;

my $LOW_QUAL_DIVISOR  = 4;
my $DEFAULT_QUAL      = 3;

#
#  Methods for reading an ACE file.
#

sub read_AS{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $num_contigs, $num_reads)=split /\s+/;
        if($id eq "AS"){
            return ($num_contigs, $num_reads);
        }
    }
    die "Could not find AS to read.\n";
}


sub read_CO{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $contig_id, $num_bases, $num_reads, $num_segments, $complementation, $sequence)=split /\s+/;

        if($id eq "CO"){
            while(<$fh>){
                chomp;
                if($_ eq ""){
                    last;
                }else{
                    $sequence.=$_;
                }
            }
            return($contig_id, $num_bases, $num_reads, $num_segments, $complementation, $sequence);
        }
    }
    die "Could not find CO to read.\n";
}


sub read_BQ{
    my $fh=shift;

    my ($id, $sequence);

    while(<$fh>){
        chomp;
        ($id)=split /\s+/;

        if($id eq "BQ"){
            while(<$fh>){
                chomp;
                if($_ eq ""){
                    last;
                }else{
                    $sequence.=$_;
                }
            }
            return($sequence);
        }
    }
    die "Could not find BQ to read.\n";
}


sub read_AF{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $read_id, $complementation, $start)=split /\s+/;
        if($id eq "AF"){
            return($read_id, $complementation, $start);
        }
    }
    die "Could not find AF to read.\n";
}


sub read_BS{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $start, $end, $read_id)=split /\s+/;
        if($id eq "BS"){
            return($start, $end, $read_id);
        }
    }
    die "Could not find BS to read.\n";
}


sub read_RD{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $read_id, $num_bases, $num_read_info_items, $num_read_tags)=split /\s+/;
        my $sequence;
        if($id eq "RD"){
            while(<$fh>){
                chomp;
                if($_ eq ""){
                    last;
                }else{
                    $sequence.=$_;
                }
            }
            return($read_id, $num_bases, $num_read_info_items, $num_read_tags, $sequence);
        }
    }
    die "Could not find RD to read.\n";
}


sub read_QA{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $qual_start, $qual_end, $clip_start, $clip_end)=split /\s+/;
        if($id eq "QA"){
            return($qual_start, $qual_end, $clip_start, $clip_end);
        }
    }
    die "Could not find QA to read.\n";
}


sub read_DS{
    my $fh=shift;
    my $id;
    while(<$fh>){
        chomp;
        my ($id)=split /\s+/;
        if($id eq "DS"){
            return("not implemented");
        }
    }
    die "Could not find DS to read.\n";
}

#
#
#

sub emitFragment ($$$$) {
    my $uid = shift;
    my $lid = shift;
    my $seq = shift;
    my $oh  = shift;

    my $len = length($seq);

    my $qvs = $seq;

    my $q = chr($DEFAULT_QUAL                   + ord("0"));
    my $l = chr($DEFAULT_QUAL/$LOW_QUAL_DIVISOR + ord("0"));

    $qvs =~ s/[^ACGT]/$l/og;
    $qvs =~ s/[ACGT]/$q/og;

    print $oh "{FRG\n";
    print $oh "act:A\n";
    print $oh "acc:$uid\n";
    print $oh "rnd:1\n";
    print $oh "sta:G\n";
    print $oh "lib:$lid\n";
    print $oh "pla:0\n";
    print $oh "loc:0\n";
    print $oh "src:\n.\n";
    print $oh "seq:\n$seq\n.\n";
    print $oh "qlt:\n$qvs\n.\n";
    print $oh "hps:\n.\n";
    print $oh "clr:0,$len\n";
    print $oh "}\n";
}

#
#
#

sub shredContig ($$$$$) {
    my $ctgId       = shift;
    my $avgCoverage = shift;
    my $sequence    = shift;
    my $libId       = shift;
    my $oh          = shift;

    my $seq_len=length($sequence);

    my @begin_shred;
    my @end_shred;

    {
        #
        #                  |*******|
        #                  |###############|
        # |-------------------------------------------------|
        #  ----------------1----------------
        #          ----------------2----------------
        #                  ----------------3----------------
        #
        #	#### represents the distance between center of read 1 and read 3
        #            [$center_range_width]
        #       **** represents the distance between centers of consective reads
        #            [$center_increments]
        #

        my $shred_len = $SHRED_READ_LENGTH;
        $shred_len = $seq_len - 50 if $seq_len < $SHRED_READ_LENGTH;

        my $num_reads=int($seq_len * $avgCoverage / $shred_len);
        my $center_range_width = $seq_len - $shred_len;

        if($num_reads==1){
            push @begin_shred, 0;
            push @end_shred, $shred_len;
        }else{
            my $center_increments = $center_range_width / ($num_reads-1);

            # Cap the number of reads we will make so that we don't get
            # redundant reads

            my $i;
            my ($prev_begin, $prev_end)=(-1,-1);
            for($i=0; $i<$num_reads; $i++){
                my $begin=$center_increments*$i;
                my $end=$begin+$shred_len;

                $begin=int($begin);
                $end=int($end);

                if($begin!=$prev_begin || $end!=$prev_end){
                    push @begin_shred, $begin;
                    push @end_shred, $end;
                    $prev_begin=$begin;
                    $prev_end=$end;
                }
            }
        }

    }

    my $num_shreds = scalar(@begin_shred);

    my $accomplished_coverage = $num_shreds * $SHRED_READ_LENGTH / $seq_len;

    # Output sequence after it has been formatted to the specified width
    my $shred_idx;
    for($shred_idx=0; $shred_idx<$num_shreds; $shred_idx++){
        my $shredded_sequence=substr($sequence,
                                     $begin_shred[$shred_idx],
                                     $end_shred[$shred_idx]-$begin_shred[$shred_idx]);

        #"/contig=$contigID\.$shred_idx " ,
        #"/target_coverage=$avgCoverage " ,
        #"/accomplished_coverage=$accomplished_coverage " ,
        #"/input_length=$seq_len " ,
        #"/range=${$begin_shred_ref}[$shred_idx]-" ,
        #       "${$end_shred_ref}[$shred_idx]\n";

        emitFragment("$libId.$ctgId.frag$shred_idx.$begin_shred[$shred_idx]-$end_shred[$shred_idx]", $libId, $shredded_sequence, $oh);
    }
}

#
#  Main
#

sub shredACE ($$) {
    my $aceFile = shift;
    my $outFile = shift;
    my $libId   = $aceFile;

    if ($aceFile =~ m/^.*\/(.*).ace/) {
        $libId = $1;
    }

    my $fh = new FileHandle "< $aceFile";
    my $oh = new FileHandle "> $outFile";

    print $oh "{VER\n";
    print $oh "ver:2\n";
    print $oh "}\n";
    print $oh "{LIB\n";
    print $oh "act:A\n";
    print $oh "acc:$libId\n";
    print $oh "ori:U\n";
    print $oh "mea:0.0\n";
    print $oh "std:0.0\n";
    print $oh "src:\n";
    print $oh ".\n";
    print $oh "nft:1\n";
    print $oh "fea:\n";
    print $oh "doNotOverlapTrim=1\n";
    print $oh ".\n";
    print $oh "}\n";

    my ($num_contigs, $num_reads)=read_AS($fh);

    my $contig_idx;
    for($contig_idx=0; $contig_idx<$num_contigs; $contig_idx++){

        my %read_position_hash;

        my ($contig_id, $num_consensus_bases, $num_reads, $num_segments, $complementation, $consensus_sequence) = read_CO($fh);

        my @coverage_array;
        my $i;

        # Initialize Coverage Array
        for($i=0; $i<$num_consensus_bases; $i++){
            $coverage_array[$i]=0;
        }

        my $quality=read_BQ($fh);

        my $read_idx;
        for($read_idx=0; $read_idx<$num_reads; $read_idx++){
            my ($read_id, $complementation, $consensus_start_pos)=read_AF($fh);
            $read_position_hash{$read_id}=$consensus_start_pos;
        }

        my ($base_line_start, $base_line_end, $base_line_read_id)=read_BS($fh);

        for($read_idx=0; $read_idx<$num_reads; $read_idx++){
            my ($read_id, $num_padded_bases, $num_read_info_items, $num_read_tags, $read_sequence)= read_RD($fh);
            my ($qual_start, $qual_end, $align_start, $align_end)=read_QA($fh);
            my $startPos = $read_position_hash{$read_id};

            my $begin = $align_start + $startPos - 1;
            my $end   = $align_end   + $startPos - 1;

            for($i=$begin; $i<$end; $i++){
                $coverage_array[$i]++;
            }
            my ($null)=read_DS($fh);
        }


        my $in_deep_enough=0;
        my @sub_contig_begin_arr;
        my @sub_contig_end_arr;

        # Keep track of where we go into deep coverage region from low coverage regions
        for($i=0; $i<$num_consensus_bases; $i++){
            if($coverage_array[$i]>$MIN_COVERAGE && !$in_deep_enough){
                push @sub_contig_begin_arr, $i;
                $in_deep_enough=1;
            }
            if($coverage_array[$i]<=$MIN_COVERAGE && $in_deep_enough){
                push @sub_contig_end_arr, ($i);
                $in_deep_enough=0;
            }
        }

        if($in_deep_enough){
            push @sub_contig_end_arr, ($i);
        }

        for($i=0; $i<=$#sub_contig_begin_arr; $i++){
            # Sum up coverage for each sub contig
            my $cov_idx;
            my $cov_sum=0;
            for($cov_idx=$sub_contig_begin_arr[$i];
                $cov_idx<$sub_contig_end_arr[$i];
                $cov_idx++){
                $cov_sum+=$coverage_array[$cov_idx];
            }

            # Compute average coverage depth

            my $sub_seq_len=$sub_contig_end_arr[$i]-$sub_contig_begin_arr[$i];
            my $avg_cov = $cov_sum / $sub_seq_len;

            if($num_reads > $MIN_READS && $sub_seq_len>=$MIN_CONTIG_SIZE){
                my $sub_contig_seq  = substr($consensus_sequence,
                                             $sub_contig_begin_arr[$i],
                                             $sub_seq_len);

                # Remove padding
                $sub_contig_seq=~s/\*//g;

                shredContig($contig_id, $avg_cov, $sub_contig_seq, $libId, $oh);
            }
        }
    }

    print $oh "{VER\n";
    print $oh "ver:1\n";
    print $oh "}\n";
}

#
#  For standalone use
#

#die "usage: $0 file.ace > file.frg\n" if (scalar(@ARGV) == 0);
#shredACE($ARGV[0], "a.frg");
#exit();
