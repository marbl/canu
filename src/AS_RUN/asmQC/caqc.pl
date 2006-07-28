#!/usr/local/bin/perl

# $Id: caqc.pl,v 1.6 2006-07-28 15:29:00 moweis Exp $
#
# This program reads a Celera .asm file and produces aggregate information
# about the assembly
#
# Written by Mihai Pop, Martin Shumway
#
#

use strict;
use IO::File;
use File::Basename;
use Statistics::Descriptive;
use TIGR::Foundation;
use TIGR::AsmLib;

use Time::HiRes;

my $MY_VERSION = " Version 2.11 (Build " . (qw/$Revision: 1.6 $/ )[1] . ")";

# Constants
my $MINQUAL = 20;
my $MINCONTIG = 10000;

my $MY_HELPTEXT = qq~
Generate quality statistics from the specified Celera assembly .asm file.  

  caqc  <prefix>.asm  [options]

    prefix.asm   The Celera .asm file
    options:
      -minqual   <n>   Minimum quality value threshhold to report as bad 
                       quality (default $MINQUAL)
      -mincontig <n>   Minimum contig size to report as a big contig
                       (default $MINCONTIG) 
      -d <c>           Specify the record delimiter tag<c>value
      -s <c>           Specify the list delimiter el[<c>el ...]
      -silent          Do not emit output to the stdout
      -g <n>           Genome size used in the calculation of N50 numbers
                       (default: TotalBasesInContigs)

caqc produces output on stdout as well as a <prefix>.qc file in the 
current directory.  Output is in a tag-value format organized in a
hierarchical manner compatible with the INI style for config files.  This
file can be read using TIGR::ConfigFile (run perldoc TIGR::ConfigFile). 
Delimiters can be controlled by using the tag-value separator specified in
the -d option, and the minor delimiter separating list values by using the 
-l option.  By default both are set to whitespace for human readability
in the stdout output, and to = and comma respectively in the .qc file.

See also:
    http://intranet.tigr.org/software_docs/CeleraAssembler.shtml
    TIGR::ConfigFile 
    ~;

my $base = new TIGR::Foundation;

if (! defined $base){
    print STDERR "Nasty error, hide!\n";
    exit(1);
}

# Global delimiters
my $DEFAULT_TAG_DELIM = '=';
my $DEFAULT_FIELD_DELIM = ',';
my $d = $DEFAULT_TAG_DELIM;   # tag-value delimiter
my $s = $DEFAULT_FIELD_DELIM;   # list separator
my $silent = 0;    


MAIN:
{
    $base->setHelpInfo($MY_HELPTEXT);
    $base->setVersionInfo($MY_VERSION);
    my $infile;
    my $delim = undef;
    my $sep = undef;
    my $genomesize;
    
    my $err = $base->TIGR_GetOptions("i=s" => \$infile,
				          "minqual=i" => \$MINQUAL,
				          "mincontig=i" => \$MINCONTIG,
                                     "d=s" => \$delim,
                                     "s=s" => \$sep,
                                     "silent!" => \$silent,
				          "g=i" => \$genomesize
				     );
    
    if ($err == 0){
	$base->bail("Command line parsing failed.  See -h option");
    }
    
    $d = (defined $delim)? $delim : $d;
    $s = (defined $sep)?   $sep   : $s;    

    if (! defined $infile){
	if ($#ARGV < 0){
	    $base->bail("Must specify an input file name.  See -h option");
	} else {
	    $infile = $ARGV[0];
	}
    }
    
    open(IN, $infile) or $base->bail("Cannot open $infile: ($!)");

    my ($prefix, $path, $suffix) = fileparse($infile, ".asm");
    
    my $record;
    
    my %lens;                # contig lengths
    #my %bads;                # contig bad qualities
    my %seqs;                # contig number of reads
    my %scaflen;             # sum of contig sizes in scaffold
    my %scafcontig;          # # contigs in scaffold
    my %inscaff;             # contigs in scaffolds
    my %adjscaflen;          # length of scaffold including inner-scaffold gaps
    my %readlens;            # sum of lengths of reads in contig
    my %utgFrags;            # frags in surrogate unitigs
    
    my $totalReadLength = 0; # Aggregate length of CLR of input reads used in assembly
    my $totalSeqs       = 0; # Total number of reads used in assembly
    
    my $ncontigs = 0;
    my $nunitigs = 0;
    my $numSingletons = 0;
    my $numVarRecords = 0;
    my $lenPlacedSurrogates = 0;
    my $numPlacedSurrogates = 0;
    my $numPlacedSurrogateFrags = 0;
    my $lenPlacedSurrogateFrags = 0;
    my $utgSeqs = 0;
    my $surrReadLen = 0;     # reads in surrogates
    my $contigReadLen = 0;   # reads in contigs
    my $degenReadLen = 0;    # reads in degenerates
    my $utgs = Statistics::Descriptive::Sparse->new();
    my %Results = ();
    $Results{ReadsWithChaffMate}     = 0;
    $Results{ReadsWithBadShortMate}  = 0;
    $Results{ReadsWithBothDegenMate} = 0;
    $Results{ReadsWithDegenMate}     = 0;
    $Results{ReadsWithDiffScafMate}  = 0;
    $Results{ReadsWithGoodMate}      = 0;
    $Results{ReadsWithBothChaffMate} = 0;
    $Results{ReadsWithBadLongMate}   = 0;
    $Results{ReadsWithNoMate}        = 0;
    $Results{ReadsWithOuttieMate}    = 0;
    $Results{ReadsWithSurrogateMate} = 0;
    $Results{ReadsWithSameOrientMate}= 0;
    $Results{ReadsWithBothSurrMate}  = 0;
    $Results{ReadsWithUnassigned}    = 0;

    my %readLen;
    my $scaffLinkWeights = 0; # sum of links connecting scaffolds

    print STDERR "Starting read ", time(), "\n";

    my $minqualtest;
    for (my $i=0; $i<$MINQUAL; $i++) {
        $minqualtest .= chr($MINQUAL + ord('0'));
    }

    my $timecco = 0;
    my $timeutg = 0;
    my $timescf = 0;
    my $timeslk = 0;
    my $timeafg = 0;

    while ($record = getCARecord(\*IN)){

        my $type;
        if ($record =~ m/^\{?(\w+)/) {
            $type = $1;
        }
	
	if ($type eq 'CCO'){
            $timecco -= time();
            my ($type, $fields, $recs) = parseCARecord($record);
	    $ncontigs++;
	    my $contiglen = $$fields{'len'};
	    my $nreads = $$fields{'npc'};
        #my $qual = $$fields{"qlt"};
	    my $id = getCAId($$fields{'acc'});
	        
	    #my @quals = split('\n', $qual);
	    #$qual = join("", @quals);
        #$qual =~ tr/\s//d;

	    my $badquals = 0;
	    #for (my $i = 0; $i < length($qual); $i++){
            #    my $thisqual = substr($qual, $i, 1);
            #    my $qualval = ord($thisqual) - ord('0');
            #    if ($qualval < $MINQUAL){
            #        $badquals++;
            #    }
	    #} # for $i < length($qual)

            #$MINQUAL += ord('0');
            #foreach my $thisqual (split '', $qual) {
            #    if (ord($thisqual) < $MINQUAL){
            #        $badquals++;
            #    }
            #}
            #$MINQUAL -= ord('0');

            #$badquals = $qual =~ s/$minqualtest/$minqualtest/go;



	    $lens{$id} = $contiglen;
	    my $ngaps = $$fields{'cns'} =~ tr/-/-/;

	    #while ($$fields{"cns"} =~ /-/g){
            #    $ngaps++;
	    #}

	    $lens{$id} -= $ngaps;
	    #$bads{$id} = $badquals;
	    $seqs{$id} = $nreads;
        foreach my $rec (@$recs){
            my ($sid, $sfs, $srecs) = parseCARecord($rec);
            if ($sid eq 'MPS'){
                #my ($l, $r) = $$sfs{'pos'};
                my $mid = $$sfs{'mid'};
                my $rlen = $readLen{$mid};#abs($r - $l);
                $readlens{$id} += $rlen;
                if ( exists $utgFrags{$mid} ) {
                    $numPlacedSurrogateFrags++;
                    $lenPlacedSurrogateFrags += $rlen;
                }
            } elsif ($sid eq 'UPS') {
                my $typ = $$sfs{'typ'};
                if ( $typ eq 'S' ) {
                    my ($ub,$ue) = split ',',$$sfs{'pos'};
                    $lenPlacedSurrogates +=  abs($ue - $ub);
                    $numPlacedSurrogates++;
                }
            } elsif ($sid eq 'VAR') {
                $numVarRecords++;
            }
        }

        $timecco += time();
	} # if $type = CCO
	
	if ($type eq 'UTG'){
            $timeutg -= time();
            my ($type, $fields, $recs) = parseCARecord($record);
	    if ($$fields{'sta'} ne 'U' && $$fields{'sta'} ne 'N'){
            $utgs->add_data($$fields{'len'});
            $nunitigs++;
            $utgSeqs += $$fields{'nfr'};
            for (my $rd = 0 ;$rd <= $#$recs; $rd++){
                my ($sid, $sfs, $srecs) = parseCARecord($$recs[$rd]);
                if ($sid eq 'MPS'){
                    my ($l, $r) = $$sfs{'pos'};
                    my $mid = $$sfs{'mid'};
                    $surrReadLen += $readLen{$mid};#abs($r - $l);
                    $utgFrags{$mid} = 1;
                }
            }
        }
            $timeutg += time();
	}
	
	if ($type eq "SCF"){
            $timescf -= time();
            my ($type, $fields, $recs) = parseCARecord($record);
	    my $id = getCAId($$fields{"acc"});
	    my %scaffcontigs = ();
	    my $adj = 0;
	    for (my $i = 0; $i <= $#$recs; $i++){
		my ($lrec, $lfield, $lrecs) = parseCARecord($$recs[$i]);
		
		$scaffcontigs{$$lfield{"ct1"}}++;
		$scaffcontigs{$$lfield{"ct2"}}++;
		$inscaff{$$lfield{"ct1"}} = 1;
		$inscaff{$$lfield{"ct2"}} = 1;
		$adj += $$lfield{"mea"};
	    }
	        
	    $scafcontig{$id} = 0;
	    $scaflen{$id} = 0;
	    while (my ($ct, $freq) = each %scaffcontigs){
		$scafcontig{$id}++;
		$scaflen{$id} += $lens{$ct};
	    }
	        
	    $adjscaflen{$id} = $scaflen{$id} + int($adj);
            $timescf += time();
	} # if $type = SCF

	if ($type eq "SLK") { # scaffold link
            $timeslk -= time();
            my ($type, $fields, $recs) = parseCARecord($record);
	    $Results{TotalScaffoldLinks}++;
	    $scaffLinkWeights += $$fields{num};
            $timeslk += time();
	}
	
	if ($type eq 'AFG') {
            $timeafg -= time();
            my ($type, $fields, $recs) = parseCARecord($record);
	    my $id = getCAId($$fields{'acc'});
	    my $clrs = getCAId($$fields{'clr'});
	    my ($clrl,$clrr) = split /,/,$clrs;
	    my $length = $clrr - $clrl;

	    $totalReadLength += $length;
	    $readLen{$id} = $length;
	    $totalSeqs++;

	    if ($$fields{'mst'} eq 'A'){
		$Results{ReadsWithChaffMate}++;
	    } elsif ($$fields{'mst'} eq 'C'){
		$Results{ReadsWithBadShortMate}++;
	    } elsif ($$fields{'mst'} eq 'D'){
		$Results{ReadsWithBothDegenMate}++;
	    } elsif ($$fields{'mst'} eq 'E'){
		$Results{ReadsWithDegenMate}++;
	    } elsif ($$fields{'mst'} eq 'F'){
		$Results{ReadsWithDiffScafMate}++;
	    } elsif ($$fields{'mst'} eq 'G'){
		$Results{ReadsWithGoodMate}++;
	    } elsif ($$fields{'mst'} eq 'H'){
		$Results{ReadsWithBothChaffMate}++;
	    } elsif ($$fields{'mst'} eq 'L'){
		$Results{ReadsWithBadLongMate}++;
        } elsif ($$fields{'mst'} eq 'N'){
		$Results{ReadsWithNoMate}++;
	    } elsif ($$fields{'mst'} eq 'O'){
		$Results{ReadsWithOuttieMate}++;
	    } elsif ($$fields{'mst'} eq 'R'){
		$Results{ReadsWithSurrogateMate}++;
	    } elsif ($$fields{'mst'} eq 'S'){
		$Results{ReadsWithSameOrientMate}++;
	    } elsif ($$fields{'mst'} eq 'U'){
		$Results{ReadsWithBothSurrMate}++;
	    } elsif ($$fields{'mst'} eq 'Z'){
		$Results{ReadsWithUnassigned}++;
	    } else {
		$Results{ReadsWithUnassigned}++;
        }

        $numSingletons++ if $$fields{'cha'} == 1;

        $timeafg += time();
	}
	
    } # while $record
    
    print STDERR "timecco = $timecco\n";
    print STDERR "timeutg = $timeutg\n";
    print STDERR "timescf = $timescf\n";
    print STDERR "timeslk = $timeslk\n";
    print STDERR "timeafg = $timeafg\n";

    print STDERR "Init ", time(), "\n";

    # Initialize the Results hash so that we never get blanks for values

    $Results{TotalScaffolds}                 = 0;
    $Results{TotalContigsInScaffolds}        = 0;
    $Results{MeanContigsPerScaffold}         = 0.0;
    $Results{MinContigsPerScaffold}          = 0;
    $Results{MaxContigsPerScaffold}          = 0;
    $Results{TotalBasesInScaffolds}          = 0;
    $Results{MeanBasesInScaffolds}           = 0.0;
    $Results{MinBasesInScaffolds}            = 0;
    $Results{MaxBasesInScaffolds}            = 0;
    $Results{N50ScaffoldBases}               = 0;
    $Results{TotalSpanOfScaffolds}           = 0;
    $Results{MeanSpanOfScaffolds}            = 0.0;
    $Results{MinScaffoldSpan}                = 0;
    $Results{MaxScaffoldSpan}                = 0;
    $Results{N50ScaffoldSpan}                = 0;
    $Results{IntraScaffoldGaps}              = 0;
    $Results{MeanSequenceGapSize}            = 0.0;
    $Results{TotalContigsInScaffolds}        = 0;
    $Results{TotalBasesInScaffolds}          = 0;
    $Results{MeanContigSize}                 = 0.0;
    $Results{MinContigSize}                  = 0;
    $Results{MaxContigSize}                  = 0;
    $Results{N50ContigBases}                 = 0;
    $Results{TotalBigContigs}                = 0;
    $Results{BigContigLength}                = 0;
    $Results{MeanBigContigSize}              = 0.0;
    $Results{MinBigContig}                   = 0;
    $Results{MaxBigContig}                   = 0;
    $Results{BigContigsPercentBases}         = 0.0;

    $Results{TotalSmallContigs}              = 0;
    $Results{SmallContigLength}              = 0;
    $Results{MeanSmallContigSize}            = 0.0;
    $Results{MinSmallContig}                 = 0;
    $Results{MaxSmallContig}                 = 0;
    $Results{SmallContigsPercentBases}       = 0.0;
    $Results{TotalDegenContigs}              = 0;
    $Results{DegenContigLength}              = 0;
    $Results{MeanDegenContigSize}            = 0.0;
    $Results{MinDegenContig}                 = 0;
    $Results{MaxDegenContig}                 = 0;
    $Results{DegenPercentBases}              = 0.0;
    $Results{TotalReads}                     = 0;
    $Results{AvgClearRange}                  = $totalReadLength / $totalSeqs;
    $Results{ReadsInContigs}                 = 0;
    $Results{BigContigReads}                 = 0;
    $Results{SmallContigReads}               = 0;
    $Results{DegenContigReads}               = 0;
    $Results{ReadsInSurrogates}              = 0;
    $Results{SingletonReads}                 = 0;
    $Results{ChaffReads}                     = $numSingletons;
    $Results{TotalNumVarRecords}             = $numVarRecords;
    $Results{AllReads}                       = 0.0;
    $Results{ContigsOnly}                    = 0.0;
    $Results{ContigsAndDegens}               = 0.0;
    $Results{TotalSurrogates}                = $nunitigs;
    $Results{SurrogateInstances}             = $numPlacedSurrogates;
    $Results{SurrogateLength}                = $utgs->sum();
    $Results{SurrogateInstanceLength}          = $lenPlacedSurrogates;
    $Results{UnPlacedSurrReadLen}     = $surrReadLen - $lenPlacedSurrogateFrags;
    $Results{PlacedSurrogateReads}           = $numPlacedSurrogateFrags;
    $Results{MeanSurrogateSize}              = $utgs->mean();
    $Results{MinSurrogateSize}               = $utgs->min();
    $Results{MaxSurrogateSize}               = $utgs->max();
    $Results{SDSurrogateSize}                = $utgs->standard_deviation();
    $Results{RangeSurrogateSize}             = $utgs->sample_range();

    if (! exists $Results{TotalScaffoldLinks}){
	$Results{TotalScaffoldLinks} = 0;
	$Results{MeanScaffoldLinkWeight} = 0;
    } else {
	$Results{MeanScaffoldLinkWeight} = 
	    $scaffLinkWeights / $Results{TotalScaffoldLinks};
    }
    
    my $totlen = 0;
    #my $totbads = 0;
    my $totseqs = 0;
    
    my $biglen = 0;
    #my $bigbads = 0;
    my $bigseqs = 0;
    my $nbigs = 0;
    
    my $maxcontig = 0;
    
    for (my $i = 0; $i < $ncontigs; $i++){
	my ($id, $len) = each(%lens);
	
	my $number;
	my $min;
	my $max;
	my $totlen;
	my $nseq;
	
	if (exists $inscaff{$id}){ # the good guys
	    if ($len >= $MINCONTIG) { # the good big guys
		$number = "TotalBigContigs";
		$min = "MinBigContig";
		$max = "MaxBigContig";
		$totlen = "BigContigLength";
		$nseq = "BigContigReads";
	    } else { # the small good guys
		$number = "TotalSmallContigs";
		$min = "MinSmallContig";
		$max = "MaxSmallContig";
		$totlen = "SmallContigLength";
		$nseq = "SmallContigReads";
	    }
	    $contigReadLen += $readlens{$id};
	} else { # the chaff
	    $number = "TotalDegenContigs";
	    $min = "MinDegenContig";
	    $max = "MaxDegenContig";
	    $totlen = "DegenContigLength";
	    $nseq = "DegenContigReads";
	    $degenReadLen += $readlens{$id};
	}

        $Results{$max} = max($len, $Results{$max});
        $Results{$min} = min($len, $Results{$min});
	$Results{$number}++;
	$Results{$totlen} += $len;
	$Results{$nseq} += $seqs{$id};
    } # for $i < $ncontigs

    $Results{'TotalContigs'} = 0 +
	$Results{'TotalBigContigs'} + $Results{'TotalSmallContigs'};
    $Results{'TotalBasesInContigs'} = 0 +
	$Results{'BigContigLength'} + $Results{'SmallContigLength'};
    $Results{'MeanContigSize'} = ($Results{'TotalContigs'} > 0)? 
	$Results{'TotalBasesInContigs'} * 1.0 
	/ $Results{'TotalContigs'} :
        0;

    $Results{'ReadsInContigs'} = 0 + 
	$Results{'BigContigReads'} + $Results{'SmallContigReads'};
    $Results{'MinContigSize'} = min($Results{'MinSmallContig'}, $Results{'MinBigContig'});
    $Results{'MaxContigSize'} = max($Results{'MaxBigContig'}, $Results{'MaxSmallContig'});
    
    $Results{'MeanBigContigSize'} = ($Results{'TotalBigContigs'} > 0)?
	$Results{'BigContigLength'} * 1.0 
	/ $Results{'TotalBigContigs'} :
        0;

    $Results{'MeanSmallContigSize'} = ($Results{'TotalSmallContigs'})? 
	$Results{'SmallContigLength'} * 1.0 
	/ $Results{'TotalSmallContigs'} :
        0;
    
    $Results{'BigContigsPercentBases'} = ($Results{'TotalBasesInContigs'} > 0)? 
	$Results{'BigContigLength'} * 100.0 
	/ $Results{'TotalBasesInContigs'} :
        0;

    $Results{'SmallContigsPercentBases'} = ($Results{'TotalBasesInContigs'} > 0)?
	$Results{'SmallContigLength'} * 100.0 
	/ $Results{'TotalBasesInContigs'} :
        0;
    
    $Results{'MeanDegenContigSize'} = ($Results{'TotalDegenContigs'})?
	$Results{'DegenContigLength'} * 1.0 
	/ $Results{'TotalDegenContigs'} :
        0;
    
    $Results{'DegenPercentBases'} = ($Results{'TotalBasesInContigs'} > 0)?
	$Results{'DegenContigLength'} * 100.0 
	/ $Results{'TotalBasesInContigs'} :
        0;

# compute N50 values and top 5 guys
    my @sortContig = sort {$lens{$b} <=> $lens{$a}} (keys %inscaff);
    my $top5contig;
    my $sum = 0;
    my $reads_tot = 0;
    my $bases_tot = 0; 
    my $top5 = @sortContig > 5 ? 5 : @sortContig;
    for (my $cc = 0; $cc < $top5; $cc++){
        $top5contig .= 
            "$cc$d$seqs{$sortContig[$cc]}$s$lens{$sortContig[$cc]}\n";
        $reads_tot += $seqs{$sortContig[$cc]};
        $bases_tot += $lens{$sortContig[$cc]};
    }
    $top5contig .= "total$d$reads_tot$s$bases_tot\n";

    my $gsz;
    if (! defined $genomesize){
        $gsz = $Results{'TotalBasesInContigs'};
    } else {
        $gsz = $genomesize;
    }
    calcNStats( 'Contig', $gsz, \%lens, \@sortContig, \%Results);
    
# do scaffold stats
    while (my ($id, $ctgs) = each %scafcontig)
    {
	$Results{'TotalScaffolds'}++;
	$Results{'TotalContigsInScaffolds'} += $scafcontig{$id};
	$Results{'IntraScaffoldGaps'} += $scafcontig{$id} - 1;

        $Results{'MinContigsPerScaffold'} = 
	    min($Results{'MinContigsPerScaffold'}, $scafcontig{$id});
        $Results{'MaxContigsPerScaffold'} = 
	    max($Results{'MaxContigsPerScaffold'}, $scafcontig{$id});
	$Results{'TotalSpanOfScaffolds'} += $adjscaflen{$id}; 

        $Results{'MinScaffoldSpan'} = 
	    min($Results{'MinScaffoldSpan'}, $adjscaflen{$id}); 
        $Results{'MaxScaffoldSpan'} = 
	    max($Results{'MaxScaffoldSpan'}, $adjscaflen{$id});

	$Results{'TotalBasesInScaffolds'} += $scaflen{$id};

        $Results{'MinBasesInScaffolds'} = 
	    min($Results{'MinBasesInScaffolds'}, $scaflen{$id});
        $Results{'MaxBasesInScaffolds'} =
	    max($Results{'MaxBasesInScaffolds'}, $scaflen{$id});
    } # while each %scafcontigs

    $Results{'MeanContigsPerScaffold'} = 
	($Results{'TotalScaffolds'} > 0)
	? $Results{'TotalContigsInScaffolds'} * 1.0 
	/ $Results{'TotalScaffolds'}
    : 0.0;
    $Results{'MeanBasesInScaffolds'} = 
	($Results{'TotalScaffolds'} > 0)
	? $Results{'TotalBasesInScaffolds'} * 1.0 
	/ $Results{'TotalScaffolds'}
    : 0.0;
    $Results{'MeanSpanOfScaffolds'} = 
	($Results{'TotalScaffolds'} > 0)
	? $Results{'TotalSpanOfScaffolds'} * 1.0 
	/ $Results{'TotalScaffolds'}      
    : 0.0;
    $Results{'MeanSequenceGapSize'} = 
	($Results{'IntraScaffoldGaps'} > 0)
	? ($Results{'TotalSpanOfScaffolds'} 
	   - $Results{'TotalBasesInScaffolds'}) * 1.0 
	       / $Results{'IntraScaffoldGaps'}
    : 0.0;
    
# compute N50 values and top 5 guys
    my @sortScaff = sort {$scaflen{$b} <=> $scaflen{$a}} (keys %scaflen);
    my $top5scaff;
    $sum = 0;
    
    my $contigs_tot = 0;
    my $size_tot = 0;
    my $span_tot = 0;
    my $avgContig_tot = 0;
    my $avgGap_tot = 0;

    $top5 = @sortScaff > 5 ? 5 : @sortScaff;
    for (my $ss = 0; $ss < $top5 ; $ss++){
        my $scf = $sortScaff[$ss];
        $top5scaff .= 
            "$ss$d$scafcontig{$scf}$s$scaflen{$scf}$s$adjscaflen{$scf}$s" .  
            sprintf("%.2f", $scaflen{$scf} * 1.0 / $scafcontig{$scf}) .  $s .
            sprintf("%.2f", ($scafcontig{$scf} - 1 > 0)?
                    (($adjscaflen{$scf} - $scaflen{$scf}) * 1.0 / ($scafcontig{$scf} - 1)) :
                    0.0
                   ) . 
            "\n"
            ;

        $contigs_tot += $scafcontig{$scf};
        $size_tot += $scaflen{$scf};
        $span_tot += $adjscaflen{$scf};
        $avgContig_tot += sprintf("%.2f", $scaflen{$scf} * 1.0 / $scafcontig{$scf});
        $avgGap_tot += sprintf("%.2f", ($scafcontig{$scf} - 1 > 0)?
                (($adjscaflen{$scf} - $scaflen{$scf}) * 1.0 / ($scafcontig{$scf} - 1)) :
                0.0);
    }
    $top5scaff .=  "total$d$contigs_tot$s$size_tot$s$span_tot$s".
                   sprintf("%.2f", $avgContig_tot). $s . 
                   sprintf("%.2f", $avgGap_tot)."\n";

	if (! defined $genomesize){
	    $gsz = $Results{'TotalBasesInScaffolds'};
	} else {
	    $gsz = $genomesize;
	}
    calcNStats( 'Scaffold', $gsz, \%scaflen, \@sortScaff, \%Results);
    
    $Results{'AllReads'} = 
	($gsz > 0)
	? sprintf("%0.2f", $totalReadLength / $gsz)
	: 0.0;

    $Results{'ContigsOnly'} = ($Results{'TotalBasesInScaffolds'} > 0) ? 
	($contigReadLen / $Results{'TotalBasesInScaffolds'}) : 
	0.0;

    $Results{'ContigAndSurrogates'} = ($Results{'TotalBasesInScaffolds'} > 0) ? 
	(($Results{'UnPlacedSurrReadLen'} + $contigReadLen) / $Results{'TotalBasesInScaffolds'}) : 
	0.0;

#    print STDERR "surrReadLen = $surrReadLen\ncontigReadLen = $contigReadLen\ndegenReadLen = $degenReadLen\ntotalReadLen = $totalReadLength\n";

    $Results{'ContigsAndDegens'} = ($Results{'TotalBasesInScaffolds'} +
				    $Results{'DegenContigLength'} > 0) 
	? (($Results{'UnPlacedSurrReadLen'} + $contigReadLen + $degenReadLen) /
	   ($Results{'TotalBasesInScaffolds'} + $Results{'DegenContigLength'}))
	: 0.0;
    
    
# Surrogates
    $Results{'ReadsInSurrogates'}       = $utgSeqs;
    
# Coverage
    $Results{'TotalReads'}              = $totalSeqs;
    $Results{'SingletonReads'}          = $Results{'TotalReads'} 
    - $Results{'BigContigReads'} 
    - $Results{'SmallContigReads'} 
    - $Results{'DegenContigReads'} 
    - $Results{'ReadsInSurrogates'}
    + $Results{'PlacedSurrogateReads'}; # Should be in contigs
    
# Parameter settings
    $Results{'MinBigContigSizeParm'}    = $MINCONTIG;
    $Results{'MinQualityParm'}          = $MINQUAL;
    
# Emit the results
    my $fh = new IO::File("> $prefix.qc") 
	or $base->bail("Could not open $prefix.qc ($!)");
    
    print "[Scaffolds]\n" if (! $silent);
    $fh->print("[Scaffolds]\n");
    printl('TotalScaffolds',          \%Results, $fh);
    printl('TotalContigsInScaffolds', \%Results, $fh);
    printlf('MeanContigsPerScaffold', \%Results, $fh);
    printl('MinContigsPerScaffold',   \%Results, $fh);
    printl('MaxContigsPerScaffold',   \%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");
    printl('TotalBasesInScaffolds',   \%Results, $fh);
    printlf('MeanBasesInScaffolds',   \%Results, $fh);
    printl('MinBasesInScaffolds',     \%Results, $fh);
    printl('MaxBasesInScaffolds',     \%Results, $fh);
    for my $val (@{$Results{'NScaffoldBases'}}) {
        printf $fh "N%2d%s=%d\n",$val->[0]*100,'ScaffoldBases',$val->[1];
    }
    for my $val (@{$Results{'IncrScaffoldBases'}}) {
        print $fh "ScaffoldAt$val->[0]=$val->[1]\n";
    }
    print "\n" if (! $silent);
    $fh->print("\n");
    printl('TotalSpanOfScaffolds',    \%Results, $fh);
    printlf('MeanSpanOfScaffolds',    \%Results, $fh);
    printl('MinScaffoldSpan',         \%Results, $fh);
    printl('MaxScaffoldSpan',         \%Results, $fh);
    printl('IntraScaffoldGaps',       \%Results, $fh);
    printlf('MeanSequenceGapSize',    \%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");
    
    print "[Top5Scaffolds${d}contigs${s}size${s}span${s}avgContig${s}avgGap)]\n" if (! $silent);
    $fh->print("[Top5Scaffolds".$d."contigs".$s."size".$s."span".$s."avgContig".$s."avgGap]\n"
	       );
    print $top5scaff if (! $silent);
    $fh->print($top5scaff);
    print "\n" if (! $silent);
    $fh->print("\n");
    
    print "[Contigs]\n" if (! $silent);
    $fh->print("[Contigs]\n");
    printl('TotalContigsInScaffolds', \%Results, $fh);
    printl('TotalBasesInScaffolds',   \%Results, $fh);
    printl('TotalNumVarRecords',      \%Results, $fh);
    printlf('MeanContigSize',         \%Results, $fh);
    printl('MinContigSize',           \%Results, $fh);
    printl('MaxContigSize',           \%Results, $fh);
    for my $val (@{$Results{'NContigBases'}}) {
        printf $fh "N%2d%s=%d\n",$val->[0]*100,'ContigBases',$val->[1];
    }
    for my $val (@{$Results{'IncrContigBases'}}) {
        print $fh "ContigAt$val->[0]=$val->[1]\n";
    }
    print "\n" if (! $silent);
    $fh->print("\n");
    
    print "[BigContigs_greater_$MINCONTIG]\n" if (! $silent);
    $fh->print("[BigContigs_greater_$MINCONTIG]\n");
    printl('TotalBigContigs',         \%Results, $fh);
    printl('BigContigLength',         \%Results, $fh);
    printlf('MeanBigContigSize',      \%Results, $fh);
    printl('MinBigContig',            \%Results, $fh);
    printl('MaxBigContig',            \%Results, $fh);
    printlf('BigContigsPercentBases', \%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");

    print "[SmallContigs]\n" if (! $silent);
    $fh->print("[SmallContigs]\n");
    printl('TotalSmallContigs',       \%Results, $fh);
    printl('SmallContigLength',       \%Results, $fh);
    printlf('MeanSmallContigSize',    \%Results, $fh);
    printl('MinSmallContig',          \%Results, $fh);
    printl('MaxSmallContig',          \%Results, $fh);
    printlf('SmallContigsPercentBases',\%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");
    
    print "[DegenContigs]\n" if (! $silent);
    $fh->print("[DegenContigs]\n");
    printl('TotalDegenContigs',       \%Results, $fh);
    printl('DegenContigLength',       \%Results, $fh);
    printlf('MeanDegenContigSize',    \%Results, $fh);
    printl('MinDegenContig',          \%Results, $fh);
    printl('MaxDegenContig',          \%Results, $fh);
    printlf('DegenPercentBases',      \%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");

    print "[Top5Contigs${d}reads${s}bases)]\n" if (! $silent);
    $fh->print("[Top5Contigs"  .$d. "reads" .$s. "bases]\n");
    print $top5contig if (! $silent);
    $fh->print($top5contig);
    print "\n" if (! $silent);
    $fh->print("\n");

    print "[Surrogates]\n" if (! $silent);
    $fh->print("[Surrogates]\n");
    printl('TotalSurrogates',         \%Results, $fh);
    printl('SurrogateInstances',      \%Results, $fh);
    printl('SurrogateLength',         \%Results, $fh);
    printl('SurrogateInstanceLength',   \%Results, $fh);
    printl('UnPlacedSurrReadLen',     \%Results, $fh);
    printl('MinSurrogateSize',        \%Results, $fh);
    printl('MaxSurrogateSize',        \%Results, $fh);
    printlf('MeanSurrogateSize',      \%Results, $fh);
    printlf('SDSurrogateSize',        \%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");

    print "[Mates]\n" if (! $silent);
    $fh->print("[Mates]\n");
    printp('ReadsWithNoMate',         \%Results, $fh, $totalSeqs);
    printp('ReadsWithGoodMate',       \%Results, $fh, $totalSeqs);
    printp('ReadsWithBadShortMate',   \%Results, $fh, $totalSeqs);
    printp('ReadsWithBadLongMate',    \%Results, $fh, $totalSeqs);
    printp('ReadsWithSameOrientMate', \%Results, $fh, $totalSeqs);
    printp('ReadsWithOuttieMate',     \%Results, $fh, $totalSeqs);
    printp('ReadsWithBothChaffMate',  \%Results, $fh, $totalSeqs);
    printp('ReadsWithChaffMate',      \%Results, $fh, $totalSeqs);
    printp('ReadsWithBothDegenMate',  \%Results, $fh, $totalSeqs);
    printp('ReadsWithDegenMate',      \%Results, $fh, $totalSeqs);
    printp('ReadsWithBothSurrMate',   \%Results, $fh, $totalSeqs);
    printp('ReadsWithSurrogateMate',  \%Results, $fh, $totalSeqs);
    printp('ReadsWithDiffScafMate',   \%Results, $fh, $totalSeqs);
    printp('ReadsWithUnassigned',     \%Results, $fh, $totalSeqs);
    printl('TotalScaffoldLinks',      \%Results, $fh);
    printlf('MeanScaffoldLinkWeight', \%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");

    print "[Reads]\n" if (! $silent);
    $fh->print("[Reads]\n");
    printl('TotalReads',              \%Results, $fh);
    printlf('AvgClearRange',           \%Results, $fh);
    printp('ReadsInContigs',          \%Results, $fh, $totalSeqs);
    printp('BigContigReads',          \%Results, $fh, $totalSeqs);
    printp('SmallContigReads',        \%Results, $fh, $totalSeqs);
    printp('DegenContigReads',        \%Results, $fh, $totalSeqs);
    printp('ReadsInSurrogates',       \%Results, $fh, $totalSeqs);
    printp('PlacedSurrogateReads',    \%Results, $fh, $totalSeqs);
    printp('SingletonReads',          \%Results, $fh, $totalSeqs);
    printp('ChaffReads',              \%Results, $fh, $totalSeqs);
    print "\n" if (! $silent);
    $fh->print("\n");

    print "[Coverage]\n" if (! $silent);
    $fh->print("[Coverage]\n");
    printlf('ContigsOnly',            \%Results, $fh);
    printlf('ContigAndSurrogates',   \%Results, $fh);
    printlf('ContigsAndDegens',       \%Results, $fh);
    printlf('AllReads',               \%Results, $fh);
    print "\n" if (! $silent);
    $fh->print("\n");
    
    $fh->close() or $base->bail("Could not close $prefix.qc ($!)"); 
    
    exit(0);
}

# Emit the results to the terminal, and if prefix specified, to the
# prefix.qc file, as a config-readable output.
#
sub printl($$$)
{
    my $tag = shift;
    my $rh_value = shift;
    my $rfh = shift;
    
    my $s1 = sprintf("%-32s", $tag);
    print(sprintf("%s%-10s", $s1, $$rh_value{$tag}) . "\n") if (! $silent); 
    $rfh->print("$tag=" . $$rh_value{$tag} . "\n");   
}

# Emit the results to the terminal, and if prefix specified, to the
# prefix.qc file, as a config-readable output.  Float form.
#
sub printlf($$$)
{
    my $tag = shift;
    my $rh_value = shift;
    my $rfh = shift;
    
    my $s1 = sprintf("%-32s", $tag);
    print(sprintf("%s%-0.2f", $s1, $$rh_value{$tag}) . "\n") if (! $silent); 
    $rfh->print("$tag=" . sprintf("%0.2f", $$rh_value{$tag}) . "\n");   
}

# output percentage of total reads

sub printp($$$$)
{
    my ($tag,$rh_value,$rfh,$tot) = @_;
    
    my $val = sprintf "%s(%.2f%%)", $$rh_value{$tag}, $$rh_value{$tag}/$tot*100;
    printf("%-32s%-10s\n", $tag, $val) if (! $silent); 
    $rfh->print("$tag=$val\n");   
}

sub max($$)
{
    my ($a,$b) = @_;
    return ($a > $b)? $a : $b;
}

sub min($$)
{
    my ($a,$b) = @_;
    return $a if ($a != 0 && $b == 0);
    return $b if ($a == 0 && $b != 0);
    return ($a < $b)? $a : $b;
}

sub calcNStats($$$$$) {
    my ($type,$gsz,$lens,$sortContig,$Results) = @_;
    my $prevSz = 0;
    my $prevIncr = 0;
    my $incrSize = 1000000;
    my @N = (.25,.5,.75);
    my $sum = 0;
    foreach my $cc (@$sortContig) {
        my $len = $lens->{$cc};
        $sum += $len;
# use while loops instead of if's since big contig/scaffold can span divisions 
        while (@N && $sum > $gsz * $N[0]) {
            push @{$Results->{"N${type}Bases"}},[shift(@N),$len];
        }
        while (int($sum / $incrSize) > $prevIncr) {
            $prevIncr++;
            if ($len != $prevSz ) {
                push @{$Results->{"Incr${type}Bases"}},[$prevIncr*$incrSize,$len];
                $prevSz = $len;
            }
        }
    }
}
