#!/usr/local/bin/perl


# $Id: ca2ace.pl,v 1.5 2007-07-31 20:04:46 eliv Exp $
#
# Converts from a Celera .asm file to a new .ace file
#

use strict;
use DB_File;
use FindBin;
use lib "$FindBin::Bin";

use TIGR::Foundation;
use TIGR::AsmLib;

my $base = new TIGR::Foundation;
my $GREP = '/bin/grep';

if (! defined $base){
    die ("Foundation cannot be created.  FATAL!\n");
}

my $VERSION = '$Revision: 1.5 $ ';
$base->setVersionInfo($VERSION);

my $HELPTEXT = q~
   ca2ace [opts] [infile]
       
   Options:   
        -i <infile>.asm TIGR .asm file.  Must end in .asm
        -o <outfile>    Output file, by default <infile>.ace
        -c <chromat_dir> Location of the chromatograms
        -p <phd_dir>    Location of the PHD directory
        -m <map>.dbm    FragId to Trace info dbm map file.
        -t <fidToTI>    FragId to Trace name or TI file.
        -s <seqPos>.dbm FragId index into .frg file dbm.
        -v contaminant scaffold id file to filter out
    ~;

$base->setHelpInfo($HELPTEXT);

my $infile;
my $outfile;
my $frgfile;
my ($fidToInfoFile, $fidToTIFile,$seqPosFile, $contamFile);

my $err = $base->TIGR_GetOptions(
				 "i=s" => \$infile,
				 "o=s" => \$outfile,
                 "m=s" => \$fidToInfoFile,
                 "s=s" => \$seqPosFile,
                 "t=s" => \$fidToTIFile,
                 "v=s" => \$contamFile
				 );

if (! $err){
    $base->bail("Error processing options!");
}

if (! defined $infile){
    $infile = $ARGV[0];
}

if (defined $infile && $infile !~ /\.*asm$/){
    $base->bail("Input file must end in .asm");
}

if (! defined $infile){
    $base->bail("You must provide an input file\n");
}
my %fidToInfo;
if (defined $fidToInfoFile) {
    tie %fidToInfo, 'DB_File', $fidToInfoFile, O_RDONLY;
}
my %fidToTraceName;
if (defined $fidToTIFile) {
    open(TI,"<$fidToTIFile") || die "Can't read $fidToTIFile";
    while(<TI>) {
        my ($fid,$tname,$centerStuff) = split;
        $fidToTraceName{ $fid } = $tname;
    }
    close TI;
}

my $prefix;

$infile =~ /(.*)\.t*asm$/;
$prefix = $1;

if (! defined $outfile){
    $outfile = "$prefix.ace";
}

$frgfile = "$prefix.frg";

my %vectors;
if ( defined $contamFile ) {
    open(VEC,"<$contamFile") || die "Can't read $contamFile";
    while(<VEC>) {
        chomp;
        $vectors{ $_ } = 1;
    }
    close VEC;
    open(VEC, "|bzip2 >$outfile.contaminant.bz2") ||
       $base->bail("Cannot open output file \"$outfile.contaminant\" : $!");
}

# we're going to pre read the DSC messages so we know we have a degen contig
# when we first read it, instead of having to wait until the end now that
# the DSC messages occur after all the CCOs
my $asmFileSize = -s $infile;
open(IN, $infile) || $base->bail("Cannot open $infile: $!\n");
# DSC messages come near the end of the .asm, might need to adjust fraction
seek IN, int( $asmFileSize * 0.95 ), 0;
my %degenCTGs;
while(<IN>) {
    if ( $_ eq "{DSC\n" ) {
        $_ = <IN>;
        die "Expected acc: got $_" unless /^acc:/;
        $_ = <IN>;
        die "Expected ctg: got $_" unless /^ctg:(\d+)/;
        $degenCTGs{ $1 } = 1;
    }
}
seek IN, 0, 0;


# Here's the process
#
# * get contig information and generate contig record in contig temp file
# * for each of the sequences generate sequence information and store coords in
# contig writing sequence records in another sequence temp file
# * generate tiling info and concatenate both temp files into an output temp
# file
#
# * when all contigs are done generate output file header and concatenate
# output temp file to output. 


# counters for the useful information
my $nContigs = 0;
my $nReads = 0;

my $fr;

open(OUT, "|bzip2 > $outfile.bz2") ||
    $base->bail("Cannot open output file \"$outfile\" : $!");

open(DEGEN, "|bzip2 >$outfile.degen.bz2") ||
    $base->bail("Cannot open output file \"$outfile.degen\" : $!");

open(SURR, "|bzip2 >$outfile.surr.bz2") ||
    $base->bail("Cannot open output file \"$outfile.surr\" : $!");

my $inAsm = 1; # are we in an assembly record
my $nseqs;     # number of sequences in this contig
my $seenseqs;  # number of sequences we've actually seen
my $contigid;  # contigId
my $seqName;
my %seqClr;    # clear range from seq file
my %seqAlnClr; # clear range from .asm file
my %seqAlnRng; # assembly range from .asm file
my %seqOff;    # sequence offset
my %rend;      # right end of sequences
my $al;        # aligned coordinates (asm_lend, rend)
my $ar;
my $end5;      # clear range
my $end3; 
my $contigLen; # # bases in the contig
my $contigSeq;

my %seqpos;
open(FRG, $frgfile) || $base->bail("Cannot open $frgfile: $!\n");
my $seekpos = tell FRG;
sub populateSeqPos($){
    my ($inFrg) = @_;

    while (my $record = getCARecord($inFrg)){
        my ($rec, $fields, $recs) = parseCARecord($record);
        if ($rec eq "FRG"){
            my $nm = $$fields{src};
            $nm =~ tr/\n/ /;
            if ($nm =~ /^\s*$/){
                $nm = $$fields{acc};
            }
            $seqpos{$$fields{acc}} = join ';',$seekpos,$nm;
        }
        $seekpos = tell FRG;
    }
}

if ( defined $seqPosFile ) {
    if (-r $seqPosFile && -s $seqPosFile) {
        tie %seqpos, 'DB_File', "$prefix.seqpos.dbm", O_RDONLY;
    } else {
        my $hashInfo = new DB_File::HASHINFO;
        #$hashInfo->{'bsize'} = 512;
        $hashInfo->{'nelem'} = 23000000;
        tie %seqpos,'DB_File',"$prefix.seqpos.dbm", O_CREAT|O_RDWR, 0644, $hashInfo;
        populateSeqPos( \*FRG );
    }
} else {
    populateSeqPos( \*FRG );
}

sub cntContigsAndReads() {
    local $/="\n{";

    my ($co,$coR);
    my ($dco,$dcoR);
    my ($sco,$scoR);
    while(<IN>) {
        if( /^UTG\nacc:\(\d+,[\s\S]+?\nsta:(\w)\n[\s\S]+?nfr:(\d+)/) {
            next unless $1 eq 'S';
            $sco++;
            $scoR += $2;
        } elsif( /^CCO\nacc:\((\d+),[\s\S]+?npc:(\d+)/) {
            if (exists $degenCTGs{ $1 }) {
                $dco++;
                $dcoR += $2;
            } else {
                $co++;
                $coR += $2;
            }
        } elsif( /^SCF\nacc:(\w+)/ ) {
            my $scf = $1;
            if (exists $vectors{ $scf }) {
                while(<IN>) {
                    if( /\nct1:(\w+)\nct2:(\w+)\n/ ) {
                        $vectors{ $1 } = 1;
                        $vectors{ $2 } = 1;
                    }
                    last if "}\n}" eq substr($_,length($_)-5,3);
                }
            }
        }
    }
    print OUT   "AS $co $coR\n\n";
    print DEGEN "AS $dco $dcoR\n\n";
    print SURR  "AS $sco $scoR\n\n";
}
# replaced by cntContigsAndReads() since we need to count 6 different things now
# little hack to allow writing of output once, takes 15s on dros but saves write/rewrite
#my $countCMD = q[perl -ne 'BEGIN{$/="\n{"} if( /^(CCO\n|UTG\n[\s\S]+?\nsta:(\w))/) {$c++;$frg=1;$c--,$frg=0 if defined $2 && $2 ne 'S'} $r++ if $frg && /^(MPS\n|UPS\ntyp:S)/; END { print "$c $r\n"}'];
#my ($numContigs,$numReads) = split ' ',`$countCMD < $infile`

cntContigsAndReads(); # write AS messages to 3 output files
seek IN,0,0;

sub printCO($$$$) {
    my ($handle, $coLine, $ctgOutRef, $seqOutRef) = @_;
    return unless @$ctgOutRef;
	print $handle $$coLine;
    print $handle @$ctgOutRef,"\n";
    print $handle @$seqOutRef,"\n";
    undef $$coLine;
    undef @$ctgOutRef;
    undef @$seqOutRef;
}
my %utgs;
my %seenutgs;
my (@ctgOut,@seqOut,$coLine);
my ($prevRecord,$prevCO) = ('','');
while (my $record = getCARecord(\*IN)){
    my ($rec, $fields, $recs) = parseCARecord($record);

    my $id = getCAId($$fields{acc});
    if ($rec eq 'AFG') {
        $seqpos{ $id } .= ";$$fields{clr}"; # assembly can change clr so use it from AFG
        next;
    }
    my $nseqs;
    $contigid = $id;

    if ( $prevRecord eq 'UTG' ) {
        printCO(\*SURR, \$coLine, \@ctgOut, \@seqOut);

    } elsif ( $prevRecord eq 'CCO' ) {
        my $handle;
        if (exists $vectors{ $prevCO } ) {
            $handle = \*VEC;
        } elsif (exists $degenCTGs{ $prevCO }) {
            $handle = \*DEGEN;
        } else {
            $handle = \*OUT;
        }
        printCO($handle, \$coLine, \@ctgOut, \@seqOut);
    }
    $prevRecord = $rec;
    $prevCO = $contigid;
         

    if ($rec eq "CCO" || $rec eq "UTG"){
        if ($rec eq "UTG" && $$fields{sta} ne "S"){
            next;
        }
	my $seq = $$fields{cns};
       $seq =~ tr/-\n/*/d;
    my $readsInTig = 0;
	if ($rec eq "UTG"){
        $readsInTig = $$fields{nfr};
	    #print STDERR "Doing unitig: $id\n";
	    $utgs{$id} = $seq;
        $utgs{$id} =~ tr/*//d; ## use ungapped

	} else {
        $readsInTig = $$fields{npc};
	    #print STDERR "Doing contig: $id\n";
	}
	$contigLen = length($seq);
    my $numGaps = $seq =~ tr/*/*/;

	my $qual = $$fields{qlt};
    $qual =~ tr/\n//d;
    my $qualLen = length($qual);
    die "$contigid seqLen $contigLen qvLen $qualLen" unless $qualLen == $contigLen;

	for (my $i = 0; $i < length($seq); $i += 50){
	    push @ctgOut, substr($seq, $i, 50), "\n";
	}
	push @ctgOut, "\n";
	$nContigs++;
	
	push @ctgOut, "BQ\n";
	my @qualvals;
	# .ace qualities are only assigned to unpadded bases
	for (my $i = 0; $i < length($qual); $i++){
	    if (substr($seq, $i, 1) ne "*"){
		push(@qualvals, ord(substr($qual, $i, 1)) - ord('0'));
	    }
	}
    die "$contigid length $contigLen gaps $numGaps numQvs ".@qualvals if
        $contigLen - $numGaps != @qualvals;
	for (my $i = 0; $i <= $#qualvals; $i+=50){
	    my $end = $i + 49;
	    if ($end > $#qualvals){
		$end = $#qualvals;
	    }
	    push @ctgOut, " ", join(" ", @qualvals[$i .. $end]), "\n";
	}
	
	push @ctgOut, "\n";

	
	%seqClr = ();
	%seqAlnClr = ();
	%seqAlnRng = ();
	%seqOff = ();
	%rend = ();

	for (my $r = 0; $r <= $#$recs; $r++){
	    my ($srec, $sfields, $srecs) = parseCARecord($$recs[$r]);
	    my $seql;
	    my $seqr;
	    my $sequence;
	    if ($srec eq "MPS" || $srec eq "UPS"){
		if ($srec eq "UPS" && $$sfields{typ} ne "S"){
		    next; # we only deal with separable unitigs
		}
		$nReads++;
		$nseqs++;
        my $fid; # fragment id
		if ($srec eq "UPS"){
            $fid = $$sfields{lid};
		    if (! exists $seenutgs{$fid}){
			$seenutgs{$fid} = "a";
		    } else {
			$seenutgs{$fid} =
			    chr(ord($seenutgs{$fid}) + 1);
		    }
		    $seqName = "surr${fid}_$seenutgs{$fid}";
		    $sequence = $utgs{$fid};
		    $seql = 0; $seqr = length($sequence);
		}else {
            $fid = $$sfields{mid};
            die "$fid not in seqpos." unless exists $seqpos{$fid};
            (my $spos,$seqName,my $sclr) = split ';',$seqpos{$fid};
		    $sequence = get_seq(\*FRG, $fid, $spos);
            $sequence =~ tr/\n//d;
		    ($seql, $seqr) = split(',', $sclr);
		}
		my @gaps = split(' ', $$sfields{del});
		my ($asml, $asmr) = split(',', $$sfields{pos});
		my $left = $seql;
		if ($asml > $asmr){
		    $sequence  = reverseComplement($sequence);
		    my $tmp = $asmr;
		    $asmr = $asml;
		    $asml = $tmp;
		    
		    $tmp = $seqr;
		    $seqr = $seql;
		    $seql = $tmp;
		    $left = length($sequence) - $seql;
		}
		# now we add gaps to the sequence
        my @outseqA;
		my $gapindex = 0;
		for (my $j = 0; $j < length($sequence); $j++){
		    my $seqj = $j - $left;# + $seql{$id} - 1; # index in untrimmed sequence
		    if ($gapindex <= $#gaps && $seqj > $gaps[$gapindex]){
			print STDERR "Weird $seqName, $seqj > $gaps[$gapindex]\n";
		    }
		    # this here is a fix for cases when the last gap index 
		    # is equal to the length of the sequence.  In this case
		    # the sequence gets an extra gap at the very end (which
		    # I might add, is completely stupid).
		    while ($gapindex <= $#gaps && $seqj == $gaps[$gapindex]){
#                       print "GS $gapindex $#gaps $seqj $gaps[$gapindex] ",
#			length($sequence), "\n";
			push @outseqA,'*';
			$gapindex++;
#                       print "GE $gapindex\n";
		    }
		    
		    push $outseqA, substr($sequence, $j, 1);
		}
		my $outseq = join '',@outseqA;

#		print "Adding seq ($seqName) = ", length($sequence), " outseq = ",
#		length($outseq), " gaps = ", $#gaps + 1, "\n";
#		print "!$sequence\n";
#		print "!$outseq\n";

		$seqAlnRng{$seqName} = sprintf("%d %d", $asml + 1, $asmr);
		$seqAlnClr{$seqName} = sprintf("%d %d", 
					       (($seql < $seqr)?$seql + 1:$seql), 
					       (($seql < $seqr)?$seqr : $seqr + 1));
		my $off = $asml;
		my $ori = ($seql > $seqr) ? "C" : "U";
		$seqOff{$seqName} = $asml;
		$rend{$seqName} = $asmr;
#		if ($ori eq "C") {$seqOff{$seqName}++;} # fix a problem in BS
		$seqOff{$seqName}++;
		if ($ori eq "C"){
		    $off -= length($sequence) - $seql;
		} else {
		    $off -= $seql;
		}

		$off++;

		push @ctgOut, "AF $seqName $ori $off\n";
		push @seqOut, "RD $seqName ", length($outseq), " 0 0\n";
		for (my $i = 0; $i <= length($outseq); $i += 50){
		    push @seqOut, substr($outseq, $i, 50), "\n";
		}
		push @seqOut, "\n";
		my $end5;
		my $end3;
		if ($ori eq "C"){
		    $end5 = length($sequence) - $seql;
		    $end3 = length($outseq) - $seqr;
		} else {
		    $end5 = $seql;
		    $end3 = length($outseq) - length($sequence) + $seqr;
		}
#		print "QA $end5 $end3 $seql $seqr ", length($sequence), " ", length($outseq), "\n";
		$end5++; #all coordinates are 1 based
		push @seqOut, sprintf("QA %d %d %d %d\n", 
				     $end5, $end3, $end5, $end3);
		my $chrmfile = "$seqName.scf";
        #my $phdfile = `/bin/echo ../phd_dir/*.$seqName.phd.1`;
		my $phdfile = "$seqName.phd.1";
		my $time;
		
		if (defined $phdfile && -r "../phd_dir/$phdfile"){
		    $time = `$GREP TIME ../phd_dir/$phdfile`;
		    $time =~ s/TIME: //;
		}
		
		if (! defined $time){
		    $time = localtime;
		} 
		
		my $dir = ($ori eq "C") ? "rev" : "forw";
		
        my $ds = "DS CHROMAT_FILE: $chrmfile PHD_FILE: $phdfile CHEM:term DYE:big TIME: $time";
        if ( exists $fidToInfo{ $fid } ) {
            my ($center,$ti,$traceName,$lib,$template,$dir) = split ',', $fidToInfo{ $fid };
            push @seqOut, "DS $traceName CENTER: $center TEMPLATE: $template DIRECTION: $dir LIBRARY: $lib\n";

        } elsif ( exists $fidToTraceName{ $fid } ) {
            push @seqOut, "$ds $fidToTraceName{$fid}\n";
        } else {
            push @seqOut, "$ds $seqName\n";
        }
    }
	} # for each subrecord
	my $prev;
	my $nBS = 0;
	foreach my $sequence ( sort {($seqOff{$a} == $seqOff{$b}) ? ($rend{$b} <=> $rend{$a}) : ($seqOff{$a} <=> $seqOff{$b})} (keys %seqOff)) {
	    if (defined $prev) {
		if ($seqOff{$sequence} - 1 < $seqOff{$prev} ||
		    $rend{$sequence} < $rend{$prev}){
		    next;
		}
		$nBS++;
		push @ctgOut, "BS $seqOff{$prev} ", $seqOff{$sequence} - 1, " $prev\n";
	    }
	    $prev = $sequence;
	}
	$nBS++;
	push @ctgOut, "BS $seqOff{$prev} $contigLen $prev\n";
	
	$coLine = "CO $contigid $contigLen $nseqs $nBS U\n";
    } # if CCO or UTG
} # while each record

# now print some info
for my $h (\*OUT, \*DEGEN, \*SURR) {
    print $h "WA{\n";
    print $h "CA_convert $VERSION\n";
    my $time = `date`;
    chomp $time;
    print $h "Run by $ENV{USER} on $time\n";
    print $h "}\n";
    close($h);
}

exit(0);

# gets a sequence by Id from a file;
sub get_seq($$$)
{
    my $file = shift;
    my $id = shift;
    my $spos = shift;

    seek $file, $spos, 0; # seek set
    my $record = getCARecord($file);
    if (! defined $record){
        print STDERR "weird error\n";
        return;
    }

    my ($rec, $fields, $recs) = parseCARecord($record);
    
    if ($rec ne "FRG"){
        print STDERR "wierd error in get_seq, expecting frg\n";
        return;
    }
    if ($$fields{acc} != $id){
        print STDERR "wierd error in get_seq, expecting $id, got $$fields{acc}\n
";
        return;
    }
    return $$fields{seq};
}
