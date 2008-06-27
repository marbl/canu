#!/usr/local/bin/perl
# (c) Copyright 2003 The Institute for Genomic Research.  All rights reserved.

###########################################################################
# $Id: ca2ta.pl,v 1.2 2008-06-27 06:29:19 brianwalenz Exp $
#
# Author: Mihai Pop / Erik Ferlanti / Abhilasha Chaudhary
#
# Description:
#
# Converts Celera Assembler output to TIGR Assembler output.
#
# Test and Debug Features:
#   Debug level 0  - All errors
#   Debug level 1  - Basic Application progress messages and warnings
#   Debug level 9  - Full Application progress messages (Verbose)
#   Debug level 10 - The works (Very verbose)
###########################################################################


# =============================== Pragmas ==================================
use strict;
use Fcntl;
use IO::Handle '_IOFBF';
use IO::File;
use File::Basename;
use TIGR::Foundation;
use TIGR::AsmLib;

# ============================  Constants ==================================
my $VERSION = "Version 3.20 (Build " . (qw/$Revision: 1.2 $/ )[1] . ")";
my $HELPTEXT = qq~
Analyze and convert Celera Assembler assembly artifacts to TIGR usable objects.

  ca2ta [options] <prefix>.asm

  options:
    -nofasta        Do not create .fasta files for contigs, placed contigs
                    and degenerate contigs output
    -nosurrfasta    Do not create .fasta files for surrogate contigs output
    -nodegenerates  Do not create .degenerates and .degenerates.fasta output
    -nofeatures     Do not create .feat output
    -o <prefix>     Output prefix (default is input prefix)
    -report         Generate a report detailing unitig information
    -unitigs        Print only the unitig records to .contig file
    -coveragestats  Run getCoverage to generate quality class and redundancy
    -nocontig       Do not create the <prefix>.contig, <prefix>.placed.contig,
                    <prefix>.degenerates.contig and <prefix>.tasm files

ca2ta converts CA output <prefix>.asm into <prefix>.fasta, <prefix>.contig, and
<prefix>.tasm similar to the corresponding outputs of TA.  The <prefix>.tasm
file is suitable for uploading to the database using aloader.  Files or links
named <prefix>.frg and <prefix>.qual must be present in the current directory.

ca2ta produces the following output files:
  <prefix>.contig             - TIGR contig(GDE) format of all contigs
  <prefix>.placed.contig      - TIGR contig(GDE) format of scaffold contigs
  <prefix>.degenerates.contig - TIGR contig(GDE) format of degenerate contigs
  <prefix>.surrogates.contig  - TIGR contig(GDE) format of surrogate contigs
  <prefix>.unitig             - Report of unitig composition of contigs
  <prefix>.unitig.config      - Config file of unitigs
  <prefix>.feat               - XML asm_feature record of unitigs
  <prefix>.tasm               - TIGR uploadable file used by timmy
  <prefix>.degenerates        - List of contigs not in a scaffold
  <prefix>.fasta              - Fasta image of all contigs
  <prefix>.placed.fasta       - Fasta image of scaffold contigs
  <prefix>.degenerates.fasta  - Fasta image of degenerate contigs
  <prefix>.surrogates.fasta   - Fasta image of surrogate contigs

~;

# Subprograms and dependent modules
# Subprograms and dependent modules
my $COVERAGE = '/usr/local/bin/getCoverage';

my @DEPENDS =
(
  "TIGR::Foundation",
  "TIGR::AsmLib",
  $COVERAGE
);

# ============================ Global Variables =============================
my %seql;                # id -> clear range correspondence
my %seqr;
my %seqpos;              # id -> file offset correspondence
my %namepos;             # name -> file offset correspondence
my %asmpos;              # position of UTG and CCO records in asm file
my %seqnames;            # id -> seqname correspondence
my %separable_unitigs;   # separable unitig ids -> occurence in nc region
my %unitigs_to_contigs;  # unitig (ids) promoted to contigs
my %tasm_subs;           # strings to substitute in .tasm file
my %unitigs;             # info on each unitig found
my %contigs;             # info on each contig
my %adjusted_contigs;    # info on adjusted contigs
my %scaffold_contigs;    # contigs appearing in a scaffold
my %feats;               # info unique type S unitigs in contigs
my %uid2iid;

my $fragfname;           # .frg file name
my $frag_fh;             # .frg file handle
my $asm_fh;              # .asm file handle
my $timesecs;
my $date;
my $who = undef;

my $tf = new TIGR::Foundation;
if (! defined $tf){
    print STDERR "the sky is falling, run away!\n";
    exit(1);
}
my $PROGRAM = $tf->getProgramInfo('name');

my $UNIQUE_FEAT_TYPE = 'CA_UNIQUE_UNITIG';
my $NONUNIQUE_FEAT_TYPE = 'CA_NONUNIQUE_UNITIG';
my $CONTIG_IDENT = 'CA_CONTIG';
my $DEGEN_IDENT = 'CA_DEGEN';
my $UNITIG_IDENT = 'CA_FREE';

my $XML_HEADER = qq~<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE FeatureSet SYSTEM "http://aserver.tigr.org/aserver/feature.dtd">
~;

my $CONTIG_SEQ_HEADER_PARSE = '^#([^[:cntrl:]\/\s\a\e]+)\(([+]?\d+)\) '.
                              '(\[RC\]|\[\]) (\d+) bases, [[:graph:]]+ '.
                              'checksum\. \{(\d+) (\d+)\} '.
                              '\<(\d+) (\d+)\>\s*$';
my $CONTIG_HEADER_PARSE = '^##([^[:cntrl:]\/\s\a\e]+) ([+-]?\d+) ([+-]?\d+) '.
                          'bases, [[:graph:]]+ checksum\.\s*$';

# ============================= MAIN ========================================
MAIN:
{
    my $contigname;
    my $orig_contigname;
    my $pcontigname;
    my $orig_pcontigname;
    my $scontigname;
    my $dcontigname;
    my $orig_dcontigname;
    my $unitigname;
    my $featurename;
    my $tigrasmname;
    my $orig_tigrasmname;
    my $tigrnuuasmname;
    my $fastafname;
    my $pfastafname;
    my $sfastafname;
    my $degenfname;
    my $degenfastafname;
    my $seek_contig_name;
    my $seek_tasm_name;
    my $seek_surr_contig_name;

    my $opref;
    my $nofasta;
    my $nosurrfasta;
    my $nodegen;
    my $nofeat;
    my $nocontig;
    my $ureport = 1;
    my $unitigs_only;
    my $coverage_stats;
    my $unitig_fh;
    my $unitig_config;

    # Set up foundation dependencies
    $tf->setHelpInfo($HELPTEXT);
    $tf->setVersionInfo($VERSION);
    $tf->addDependInfo(@DEPENDS);

    my $err = $tf->TIGR_GetOptions("nofasta"       => \$nofasta,
                                   "nosurrfasta"   => \$nosurrfasta,
                                   "nodegenerates" => \$nodegen,
                                   "nofeatures"    => \$nofeat,
                                   "nocontig"      => \$nocontig,
                                   "report!"       => \$ureport,
                                   "o=s"           => \$opref,
                                   "unitigs"       => \$unitigs_only,
                                   "coveragestats" => \$coverage_stats);
    if ($err == 0){
        $tf->bail("Command line parsing failed.  See -h option");
    }

    ### try to read and open the input .asm file
    if (isReadableFile($ARGV[0])) {
        $asm_fh = new IO::File "$ARGV[0]" or
          $tf->bail("Cannot open \"$ARGV[0]\": $!");
    } else {
        $tf->bail("You must specify a .asm file!");
    }

    $ARGV[0] =~ /(\S+)\.asm/;
    if (! defined $1){
        $tf->bail("\"$ARGV[0]\" not recognized, must have extension .asm");
    }

    if (!defined $opref){
        $opref = $1;
    }
    $contigname = "$opref.contig";
    $orig_contigname = "$opref.orig.contig";
    $pcontigname = "$opref.placed.contig";
    $orig_pcontigname = "$opref.orig.placed.contig";
    $scontigname = "$opref.surrogates.contig";
    $dcontigname = "$opref.degenerates.contig";
    $orig_dcontigname = "$opref.orig.degenerates.contig";
    $unitigname = "$opref.unitig";
    $featurename = "$opref.feat";
    $tigrasmname = "$opref.tasm";
    $orig_tigrasmname = "$opref.orig.tasm";
    $tigrnuuasmname = "$opref.nuu.tasm";
    $fragfname = "$1.frg";
    $fastafname = "$opref.fasta";
    $pfastafname = "$opref.placed.fasta";
    $sfastafname = "$opref.surrogates.fasta";
    $degenfname = "$opref.degenerates";
    $degenfastafname = "$opref.degenerates.fasta";
    $seek_contig_name = "$opref.contig.seek";
    $seek_tasm_name = "$opref.tasm.seek";
    $seek_surr_contig_name = "$opref.contig.surrogates.seek";
    $unitig_config = "$opref.unitig.config";

    $frag_fh = new IO::File "$fragfname" or
        $tf->bail("Cannot open $fragfname: $!");

    $tf->logLocal("Getting date from file &fragfname ...", 1);

    my $record;
    while ( $record = getCARecord( $frag_fh ) ) {
      my $type;
      if ( $record =~ m/^\{?(\w+)/ ) {
    	$type = $1;
      }

      my ( $type, $fields, $recs ) = parseCARecord($record);
      if ( $type eq 'BAT' ) {  #BAT is the first message of a frg
	    $timesecs = $$fields{'crt'};
      } elsif ( $type eq 'ADT' ) {  #ADT is the second message of a frg
        my ( $lrec, $lfield, $lrecs ) = parseCARecord( $$recs[0] ); #Read the ADL
        $who      = $$lfield{'who'};
      } else {
        last;
      }
    }

    $who = $ENV{USER} if ( !defined $who or $who eq '');

    my ($seconds, $minutes, $hours, $day_of_month, $month, $year,
        $wday, $yday, $isdst) = localtime($timesecs);

    my $time_type = undef;
    if($hours >= 12) {
       $time_type = "PM";
       if($hours > 12) {
          $hours = $hours - 12;
       }
    }
    else {
       $time_type = "AM";
    }
    $year = $year + 1900;
    $year =~ /\d\d(\d\d)/;
    $year = $1;
    $date = sprintf("%02d/%02d/%02d %02d:%02d:%02d",
               $month+1, $day_of_month, $year, $hours, $minutes, $seconds);
    $date = $date." $time_type";

    $frag_fh->seek(0,SEEK_SET);

    # parse the asm file and gather the relevant information
    $tf->logLocal("Parsing asm file, $ARGV[0]...", 1);
    my $record = "";
    while ( defined(my $seekpos = tell $asm_fh) and
            ($record = getCARecord($asm_fh)) ) {
        my ($rec, $fields, $recs) = parseCARecord($record);

        if ($rec eq "AFG"){ # augmented fragment
            processAugmentedFragment($fields,$recs);
        } elsif ($rec eq "UTG") { # unitig record
            processUnitig($fields,$recs,$seekpos);
        } elsif ($rec eq "CCO") { # contig record
            processContig($fields,$recs,$seekpos,$unitig_fh);
        } elsif ($rec eq "SCF") { # scaffold record
            processScaffold($fields,$recs);
        } else {
            $tf->logLocal("$rec ($$fields{acc}) is not currently parsed and" .
                            " has " . ($#$recs + 1) . " sub-records", 9);
        }
    }
    $tf->logLocal("Indexing the .frg file", 3);
    ### index the frg file
    indexFrgFile($frag_fh);

    if ($unitigs_only) {  # print only unitigs to contig file
        $tf->logLocal("Printing only unitigs to the contig file, " .
                      "$contigname", 1);
        my $contig_fh = new IO::File "> $contigname" or
           $tf->bail("Cannot create $contigname: $!");
        foreach my $unitig_id (sort keys %unitigs) {
            $tf->logLocal("Printing unitig id '$unitig_id' to " .
                          "$contigname...", 9);
            printUnitig($contig_fh, $unitig_id, 'contig');
	}
        $contig_fh->close();
    } else {
        ### process contigs and build information on separable unitigs
        $tf->logLocal("Building information on separable (Type S) unitigs", 1);
        foreach my $contig_id (sort keys %contigs) {
            checkForSeparableUnitigs($contig_id);
        }

        if ($ureport) {
            generateUnitigReport($unitigname, $unitig_config);
        }

        ### print contigs to .contig and .fasta files
        $tf->logLocal("Printing contigs to $contigname, $pcontigname," .
                      " and $dcontigname...", 1);

        my($contig_fh,$pcontig_fh,$dcontig_fh,$seek_contig_fh,
           $seek_tasm_fh,$seek_surr_contig_fh,$orig_contig_fh,
           $scontig_fh);

        if (!$nocontig ) {
           $contig_fh = new IO::File "> $contigname" or
              $tf->bail("Cannot create $contigname: $!");
           $pcontig_fh = new IO::File "> $pcontigname" or
              $tf->bail("Cannot create $pcontigname: $!");
           $dcontig_fh = new IO::File "> $dcontigname" or
              $tf->bail("Cannot create $dcontigname: $!");
        }

        $orig_contig_fh = new IO::File ">$orig_contigname" or
           $tf->bail("Cannot open \"$orig_contigname\": $!\n");
        $scontig_fh = new IO::File "> $scontigname" or
           $tf->bail("Cannot create $scontigname: $!");

        if($nocontig) {
           $seek_contig_fh = new IO::File "> $seek_contig_name" or
              $tf->bail("Cannot create file $seek_contig_name: $!");
           $seek_tasm_fh = new IO::File "> $seek_tasm_name" or
              $tf->bail("Cannot create file $seek_tasm_name: $!");
           $seek_surr_contig_fh = new IO::File "> $seek_surr_contig_name" or
              $tf->bail("Cannot create file $seek_surr_contig_name: $!");
        }

        my($fasta_fh,$sfasta_fh,$pfasta_fh,$dfasta_fh);
        if(!$nofasta ){
           $tf->logLocal("Printing all contigs to $fastafname", 1);
           $fasta_fh = new IO::File "> $fastafname" or
              $tf->bail("Cannot create $fastafname: $!");
           $tf->logLocal("Printing placed contigs to $pfastafname", 1);
           $pfasta_fh = new IO::File "> $pfastafname" or
              $tf->bail("Cannot create $pfastafname: $!");

           if(! $nodegen) {
              $tf->logLocal("Printing degenerate contigs to $degenfastafname",
                            1);
              $dfasta_fh = new IO::File "> $degenfastafname" or
                 $tf->bail("Cannot create $degenfastafname: $!");
	   }
        }

        if(!$nosurrfasta) {
           $tf->logLocal("Printing surrogate contigs to $sfastafname", 1);
           $sfasta_fh = new IO::File "> $sfastafname" or
              $tf->bail("Cannot create $sfastafname: $!");
        }

        my $prev_pos = 0;
        foreach my $contig_id (sort keys %contigs) {
           $tf->logLocal("Printing contig id '$contig_id' to " .
                         "$contigname and $fastafname...", 9);
           printContig($contig_fh, $contig_id, 'contig') if (!$nocontig);
           printContig($orig_contig_fh, $contig_id, 'contig');

           if($nocontig) {
	      my $pos = tell($orig_contig_fh);
	      $seek_contig_fh->print("$contig_id $prev_pos $pos ");
              $prev_pos = $pos;
           }

           if(!$nofasta) {
	      printFasta($fasta_fh, $contig_id);
           }

           if(exists $scaffold_contigs{$contig_id}) {
              $tf->logLocal("Printing contig id '$contig_id' to " .
                      "$pcontigname and $pfastafname...", 9);
              if($nocontig) {
	         $seek_contig_fh->print("1\n");
              }
              else {
		 printContig($pcontig_fh, $contig_id, 'contig');
              }
              printFasta($pfasta_fh, $contig_id) if (!$nofasta);
           }
           else {
              $tf->logLocal("Printing contig id '$contig_id' to " .
                            "$dcontigname and $degenfastafname...", 9);

              if($nocontig) {
		 $seek_contig_fh->print("0\n");
              }
              else {
		 printContig($dcontig_fh, $contig_id, 'contig');
              }
              printFasta($dfasta_fh, $contig_id) if ((!$nofasta) &&
                                                     (!$nodegen));
           }
        }

        my $prev_pos1 = 0;
        foreach my $contig_id (sort keys %separable_unitigs) {
           $tf->logLocal("Printing promoted unitig id '$contig_id' to " .
                         "$scontigname...", 9);
           printUnitig($scontig_fh, $contig_id, 'contig');

           if($nocontig) {
	      my $pos = tell($scontig_fh);
	      $seek_surr_contig_fh->print("$contig_id $prev_pos1 $pos\n");
              $prev_pos1 = $pos;
           }

           if(!$nosurrfasta) {
              $tf->logLocal("Printing promoted unitig id '$contig_id' to " .
                            "$sfastafname...", 9);
              printUnitigToFasta($sfasta_fh, $contig_id);
           }
        }

        if(defined $contig_fh) {
           $contig_fh->close();
        }
        if(defined $pcontig_fh) {
           $pcontig_fh->close();
        }
        if(defined $dcontig_fh) {
           $dcontig_fh->close();
        }
        if(defined $seek_contig_fh) {
           $seek_contig_fh->close();
        }
        if(defined $seek_surr_contig_fh) {
           $seek_surr_contig_fh->close();
        }
        if(defined $orig_contig_fh) {
           $orig_contig_fh->close();
        }
        if(defined $scontig_fh) {
           $scontig_fh->close();
        }
        if(defined $fasta_fh) {
           $fasta_fh->close();
        }
        if(defined $pfasta_fh) {
           $pfasta_fh->close();
        }
        if(defined $dfasta_fh) {
           $dfasta_fh->close();
	}
        if(defined $sfasta_fh) {
           $sfasta_fh->close();
	}

        if (defined $coverage_stats) {
            ### process contigs using getCoverage to generate quality class
            ### and redundancy information
            $tf->logLocal("Building coverage information on contigs...", 1);
            if(! $nocontig) {
               calculateContigInfo($contigname);
               calculateContigInfo($scontigname);
	    }
            else {
               calculateContigInfo($orig_contigname);
               calculateContigInfo($scontigname);
            }
        }

        ### print contigs to .tasm file
        $tf->logLocal("Printing contigs to $tigrasmname", 1);
        my $tasm_fh = undef;
        my $orig_tasm_fh = undef;

        if(! $nocontig) {
           $tasm_fh = new IO::File "> $tigrasmname" or
              $tf->bail("Cannot create $tigrasmname: $!");
        }
        else {
           $orig_tasm_fh = new IO::File "> $orig_tigrasmname" or
              $tf->bail("Cannot create $orig_tigrasmname: $!");
        }

        my $nuu_tasm_fh = new IO::File "> $tigrnuuasmname" or
              $tf->bail("Cannot create $tigrnuuasmname: $!");

        my $first = 1;
        my $prev_pos = 0;

        foreach my $contig_id (sort keys %contigs) {
           $tf->logLocal("Printing contig id '$contig_id' to " .
                          "$tigrasmname...", 9);
           if ($first) {
              if(! $nocontig) {
	         printContig($tasm_fh, $contig_id, 'asm', 'first');
	      }
              else {
                 printContig($orig_tasm_fh, $contig_id, 'asm', 'first');
	      }
              $first = 0;
           }
           else {
              if(! $nocontig) {
                 printContig($tasm_fh, $contig_id, 'asm');
	      }
              else {
                 printContig($orig_tasm_fh, $contig_id, 'asm');
	      }
           }

           if($nocontig) {
	      my $pos = tell($orig_tasm_fh);
	      $seek_tasm_fh->print("$contig_id $prev_pos $pos\n");
              $prev_pos = $pos;
           }
        }

        ## print out the promoted unitigs to the .tasm file
        $tf->logLocal("Printing non-unique unitigs to $tigrasmname...", 1);

        foreach my $contig_id (sort keys %unitigs_to_contigs) {
           if(! $nocontig) {
              $tf->logLocal("Printing contig id '$contig_id' to " .
                          "$tigrasmname...", 9);
              printUnitig($tasm_fh, $contig_id, 'asm');
           }
           printUnitig($nuu_tasm_fh, $contig_id, 'asm');
        }

        if(defined $tasm_fh) {
           $tasm_fh->close();
        }
        if(defined $orig_tasm_fh) {
           $orig_tasm_fh->close();
        }
        if(defined $nuu_tasm_fh) {
           $nuu_tasm_fh->close();
        }

        ### output feature and degenerates files
        if ( ! $nofeat ) {
           outputFeatures($featurename);
        }

        if ( ! $nodegen ) {
           outputDegenerates($degenfname);
        }
     }

    if(defined $frag_fh) {
       $frag_fh->close();
    }
    if(defined $asm_fh) {
       $asm_fh->close();
    }

    exit(0);
}

######################################################################
# This method checks for separable (Type S) unitigs in a contig. If
# any occur, the contig is added to the adjusted_contigs hash and
# deleted from the contigs hash. It accepts a contig id and returns
# nothing.
######################################################################
sub checkForSeparableUnitigs {
    my $contig_id = shift;

    $tf->logLocal("Processing contig id '$contig_id' for " .
                  "info on separable unitigs...", 4);
    my $separables = 0;
    my $nseqs = 0;
    my $ra_utgs = $contigs{$contig_id}->{'utgs'};
    my $contig_len = $contigs{$contig_id}->{'len'};
    foreach my $utgr (@$ra_utgs) {
        my($utg,$upos) = split /~~/, $utgr;
        my $u_nseqs = $unitigs{$utg}->{'nseq'};
        my $utype = $unitigs{$utg}->{'type'};

        if ($utype eq 'S' and $separable_unitigs{$utg} > 1) {
            $tf->logLocal("Found non-unique type S unitig $utg in " .
                          "contig $contig_id with $u_nseqs seqs", 4);
            $separables++;
            $unitigs_to_contigs{$utg}++;
            $u_nseqs = 0;
        } elsif ($utype eq 'S') {
            $tf->logLocal("Found unique type S unitig $utg in " .
                          "contig $contig_id with $u_nseqs seqs", 4);
            $separables++;
        } else {
            $tf->logLocal("Found type $utype unitig $utg in " .
                          "contig $contig_id with $u_nseqs seqs", 4);
        }
        $nseqs += $u_nseqs;
        $tf->logLocal("Adding $u_nseqs seqs to contig. Total = $nseqs", 4);
    }

    if ($separables) {
        $adjusted_contigs{$contig_id}->{'utgs'} = $ra_utgs;
        $adjusted_contigs{$contig_id}->{'nseq'} = $nseqs;
        $adjusted_contigs{$contig_id}->{'len'} = $contig_len;
    }
}

######################################################################
# This method indexes the fragment file. It accepts as input a
# filehandle and populates a hash of file positions indexed by the
# fragment id.
######################################################################
sub indexFrgFile {
    my $fh = shift;
    $tf->logLocal("Indexing the frg file, $fragfname...", 1);
    my $seekpos = tell $fh;
    while (my $record = getCARecord($fh)) {
        my ($rec, $fields, $recs) = parseCARecord($record);
        if ($rec eq "FRG") {
            $seqpos{$$fields{acc}} = $seekpos;
            my $nm = $$fields{src};
            my @lines = split('\n', $nm);
            $nm = join(' ', @lines);
            if ($nm eq "" || $nm =~ /^\s*$/) {
                $nm = "$$fields{acc}";
            }
            $seqnames{$$fields{acc}} = $nm;
            $namepos{$nm} = $seekpos;
        }
        $seekpos = tell $fh;
    }
}

######################################################################
# This method processes a scaffold record (SCF). It accepts a
# reference to a hash containing the fields of the record and a
# reference to an array containing any sub-records. A hash containing
# the contigs in a scaffold is populated.
######################################################################
sub processScaffold {
    my $fields = shift;
    my $recs = shift;
    my $id = getCAId($$fields{'acc'});
    my $num_ctgs = $$fields{'noc'} + 1;

    $tf->logLocal("Processing scaffold $id with $num_ctgs contigs...", 4);

    # here we parse the individual unitigs aligned to the contig
    for (my $i = 0; $i <= $#$recs; $i++) {
        my($lrec,$lfield,$lrecs) = parseCARecord($$recs[$i]);
        my $contig1 = $$lfield{'ct1'};
        $tf->logLocal("Found contig $contig1 in scaffold $id", 4);
        $scaffold_contigs{$contig1}++;
        my $contig2 = $$lfield{'ct2'};
        $tf->logLocal("Found contig $contig2 in scaffold $id", 4);
        $scaffold_contigs{$contig2}++;
    }
}

######################################################################
# This method processes an augmented fragment record (AFG). It accepts
# a reference to a hash containing the fields of the record and a
# reference to an array containing any sub-records. A hash containing
# the clear left and clear right values is populated.
######################################################################
sub processAugmentedFragment {
    my $fields = shift;
    my $recs = shift;
    my $id = getCAId($$fields{acc});

    $tf->logLocal("Processing AFG $id...", 9);

    if ($#$recs != -1) {
        $tf->logLocal("Fragment $id matches " . $#$recs + 1 . " screens", 1);
    }

    # convert the external ids into internal ids
    $$fields{acc} =~ /\((\d+),(\d+)\)/;
    if (!defined $1){
        $tf->bail("Can't understand accession \"$$fields{acc}\"");
    }
    $uid2iid{$1} = $2;

    # get the new clear range
    $$fields{clr} =~ /(\d+),(\d+)/;
    if (!defined $1){
        $tf->logError("Weird clear range for fragment, $id");
    } else {
        $seql{$id} = $1;
        $seqr{$id} = $2;
    }
}

######################################################################
# This method prints a unitig record to a fasta file. It accepts a
# filehandle and a unitig id.
######################################################################
sub printUnitigToFasta {
    my $fh = shift;
    my $unitig_id = shift;

    my($rec,$fields,$recs) = getAsmRecord($unitig_id);
    my $uid = getCAId($$fields{acc});

    my $len = $$fields{len};
    my $lseq = $$fields{cns};
    my $nseq = $$fields{nfr};

    my @fields = split('\n', $lseq);
    $lseq = join('', @fields);

    print_consensus($fh, $uid, $len, $nseq, $lseq, 'fasta');
}

######################################################################
# This method prints a contig consensus to a fasta file. It accepts
# a filehandle and a contig id.
######################################################################
sub printFasta {
    my $fh = shift;
    my $contig_id = shift;
    my $nseq = shift || undef;

    my($rec,$fields,$recs) = getAsmRecord($contig_id);
    my $id = getCAId($$fields{acc});

    my $len = $$fields{len};
    my $lseq = $$fields{cns};
    if (!defined $nseq) {
        $nseq = $$fields{npc};
    }

    my @seqlines = split('\n', $lseq);
    $lseq = join('', @seqlines);

    print_consensus($fh, $id, $len, $nseq, $lseq, "fasta");
}

######################################################################
# This method processes a contig record (CCO). It accepts a reference
# to a hash containing the fields of the record, a reference to an
# array containing any sub-records, and a file position of this record
# in the .asm file. A contig hash is populated with all the relevant
# information.
######################################################################
sub processContig {
    my $fields = shift;
    my $recs = shift;
    my $seekpos = shift;
    my $unitig_fh = shift;

    my $id = getCAId($$fields{acc});

    $tf->logLocal("Gathering info for contig $id...", 9);
    $asmpos{$id} = $seekpos;

    $contigs{$id}->{'utgs'} = [];

    my $len = $$fields{len};
    my $nseq = $$fields{npc};
    $contigs{$id}->{'len'} = $len;
    $contigs{$id}->{'nseq'} = $nseq;

    # here we parse the individual unitigs aligned to the contig
    for (my $i = 0; $i <= $#$recs; $i++){
        my ($sid, $sfs, $srecs) = parseCARecord($$recs[$i]);
        if ($sid eq "UPS") {
            my $uid = $$sfs{'lid'};
            my $utype = $unitigs{$uid}->{'type'};
            my $ulen = $unitigs{$uid}->{'len'};
            if ($utype eq 'S') {
               $separable_unitigs{$uid}++;
            }
            my $u_pos = $uid . "~~" . $$sfs{'pos'};
            my $u_gaps = $$sfs{'del'};
            $tf->logLocal("Got the following gaps for unitig $uid " .
              "in contig $id: $u_gaps", 9);
            my $c_gaps = "$id~~". $$sfs{'pos'} . "~~$u_gaps";
            $unitigs{$uid}->{'gaps'} =  $u_gaps;
            push(@{$unitigs{$uid}->{'in_contig'}}, $c_gaps);
            push(@{$contigs{$id}->{'utgs'}}, $u_pos);
        }
    }
}

######################################################################
# This method prints a contig to the .contig or .tasm file. It accepts
# a filehandle, a contig id, a type (contig|asm) which describes the
# format to print, and a flag (0|1) which signifies whether this is
# the first record to print.
######################################################################
sub printContig {
    my $fh = shift;
    my $contig_id = shift;
    my $ftype = shift;
    my $first = shift or 0;

    my($rec,$fields,$recs) = getAsmRecord($contig_id);
    my $id = getCAId($$fields{acc});

    my $len = $$fields{len};
    my $lseq = $$fields{cns};
    my $nseq = $$fields{npc};

    my $seq_len = length($lseq);

    # if the length of the contig consensus is  zero do not print the contig
    if(length($lseq) == 0) {
        $tf->logLocal("WARNING: The contig $id consensus length is zero", 4);
        return;
    }
    my @fields = split('\n', $lseq);
    $lseq = join('', @fields);

    my $ident = (exists $scaffold_contigs{$contig_id}) ?
                $CONTIG_IDENT :
                $DEGEN_IDENT;

    if ($ftype eq 'asm' and $first) {
        print_consensus($fh,$id,$len,$nseq,$lseq,$ftype,$ident,'first');
    } elsif ($ftype eq 'asm') {
        print_consensus($fh,$id,$len,$nseq,$lseq,$ftype,$ident);
    } elsif ($ftype eq 'contig') {
        print_consensus($fh,$id,$len,$nseq,$lseq,$ftype);
    }

    my $ra_offsets = buildOffsets($lseq);

    $tf->logLocal("Printing aligned sequences for contig $id...", 9);
    my $tmpfile = "tmpfile";
    if(-e $tmpfile) {
       if(!unlink $tmpfile) {
          $tf->bail("Could not remove $tmpfile\n");
       }
    }
    my $tmp_fh = new IO::File ">$tmpfile" or
           $tf->bail("Cannot create $tmpfile: $!");
    my %seqs = ();
    # here we parse the individual sequences aligned to the contig
    for (my $i = 0; $i <= $#$recs; $i++){
        my ($sid, $sfs, $srecs) = parseCARecord($$recs[$i]);
        if ($sid eq "MPS" and $$sfs{typ} eq "R"){
            my ($id, $sequence, $asml, $asmr, $rc, $seqleft, $seqright) =
                processMPS($sfs);
	    $tf->logLocal("Printing sequence $seqnames{$id}.", 9);
            my $start = tell($tmp_fh);
            print_aligned($tmp_fh, $seqnames{$id}, $sequence, $asml, $rc,
                          $seqleft, $seqright,$$ra_offsets[$asml],
                          $$ra_offsets[$asmr - 1], $ftype);
            my $end = tell($tmp_fh);
            $seqs{$id}->{"offset"} = $asml;
            $seqs{$id}->{"start"} = $start;
            $seqs{$id}->{"end"} = $end;
            $seqs{$id}->{"len"} = length $sequence;

        } else {
            $tf->logLocal("Record $sid (type $$sfs{typ}) not a contig" .
                         " sequence, ignoring", 9);
        }
    }
    close $tmp_fh;
    $tmp_fh = new IO::File "$tmpfile" or
           $tf->bail("Cannot create $tmpfile: $!");
    foreach my $seq_id ( sort {
	if(($seqs{$a}->{"offset"}) < ($seqs{$b}->{"offset"})) {
	    return -1;
        }
        elsif(($seqs{$a}->{"offset"}) == ($seqs{$b}->{"offset"})) {
            if(($seqs{$a}->{"len"}) < ($seqs{$b}->{"len"})) {
	       return -1;
            }
	    elsif(($seqs{$a}->{"len"}) == ($seqs{$b}->{"len"})) {
		return 0;
	    }
            elsif(($seqs{$a}->{"len"}) > ($seqs{$b}->{"len"})) {
		return 1;
            }
        }
        elsif(($seqs{$a}->{"offset"}) > ($seqs{$b}->{"offset"})) {
	    return 1;
        } } (keys %seqs)) {

        my $start = $seqs{$seq_id}->{"start"};
        my $end = $seqs{$seq_id}->{"end"};
        my $seq_record = undef;

        if ( ! seek $tmp_fh, $start, 0) {
           $tf->bail("Trouble seeking to position $start in tmpfile");
        }

        if(!defined read($tmp_fh, $seq_record, ($end - $start))) {
           $tf->bail("Could not read the sequence information from the tmpfile");
        }
        print $fh $seq_record;
    }
    close $tmp_fh;
    if(-e $tmpfile) {
       if(!unlink $tmpfile) {
          $tf->bail("Could not remove $tmpfile\n");
       }
    }
}

########################################################################
# This method build an offset array. It accepts a sequence.
# Here the offsets array will contain one entry for each position
# in the consensus ($lseq).  That entry is the number of non-gap
# characters in the consensus, up to the current coordinate.
# These numbers follow the asm_lend, asm_rend convention and start at 1.
# Returns a reference to the offset array.
########################################################################
sub buildOffsets {
    my $lseq = shift;
    my @offsets = ();

    $#offsets = length($lseq) - 1;
    my $coord = 0;
    for (my $i = 0; $i < length($lseq); $i++){
        if (substr($lseq, $i, 1) ne "-"){
            $coord++;
        }
        $offsets[$i] = $coord;
    }
    return \@offsets;
}

######################################################################
# This method gets record from .asm file. It accepts an id which is
# a key into the asm position hash. Returns the record type, a
# reference to a hash of record fields, and a reference to an array
# of sub-records.
######################################################################
sub getAsmRecord {
    my $id = shift;

    $tf->logLocal("Grabbing record for contig/unitig $id at " .
                  "position $asmpos{$id} from $ARGV[0]", 9);
    if ( ! seek $asm_fh, $asmpos{$id}, 0) {
        $tf->bail("Trouble seeking to position $asmpos{$id}" .
                    " in $ARGV[0]");
    }
    my $record = getCARecord($asm_fh);
    if (!defined $record) {
        $tf->bail("Trouble getting record from $ARGV[0]");
    }
    my($rec,$fields,$recs) = parseCARecord($record);
    return ($rec,$fields,$recs);
}

######################################################################
# This method processes a unitig record (UTG). It accepts a reference
# to a hash of containing the fields of the record, a reference to
# an array containing the sub-records, and a file position of this
# record in the .asm file. A unitig hash is populated with the
# relevant information.
######################################################################
sub processUnitig {
    my $fields = shift;
    my $recs = shift;
    my $seekpos = shift;
    my $uid = getCAId($$fields{acc});
    $tf->logLocal("Gathering info for unitig $uid...", 9);

    $asmpos{$uid} = $seekpos;

    $unitigs{$uid}->{'len'} = $$fields{'len'};
    $unitigs{$uid}->{'type'} = $$fields{'sta'};
    $unitigs{$uid}->{'nseq'} = $$fields{'nfr'};
    $unitigs{$uid}->{'in_contig'} = [];

    if ($$fields{"sta"} eq "S"){
        $separable_unitigs{$uid} = 0;
    }

    my @seqs = ();
    # here we parse the individual sequences aligned to the unitig
    $tf->logLocal("Building list of sequences aligned to unitig $uid...", 9);
    for (my $i = 0; $i <= $#$recs; $i++){
        my ($sid, $sfs, $srecs) = parseCARecord($$recs[$i]);
        if ($sid eq "MPS") {
            my $seq_id = $$sfs{'mid'};
            $tf->logLocal("Adding seq $seq_id", 9);
            push(@seqs, $seq_id);
        }
    }
    $unitigs{$uid}->{'seqs'} = [@seqs];
}

######################################################################
# This method prints a unitig to a .contig or .tasm file. It accepts
# a filehandle, a unitig id, a type (contig|asm) which describes the
# format to print, and a flag (0|1) which tells which direction to
# print the unitig.
######################################################################
sub printUnitig {
    my $fh = shift;
    my $unitig_id = shift;
    my $ftype = shift;
    my $reverse = shift || undef;

    my($rec,$fields,$recs) = getAsmRecord($unitig_id);
    my $uid = getCAId($$fields{acc});

    my $len = $$fields{len};
    my $lseq = $$fields{cns};
    my $nseq = $$fields{nfr};

    my @fields = split('\n', $lseq);
    $lseq = join('', @fields);

    if ($reverse) {
       $lseq = reverseComplement($lseq);
    }

    if ($ftype eq 'asm') {
       print_consensus($fh, $uid, $len, $nseq, $lseq, $ftype, $UNITIG_IDENT);
    }
    elsif ($ftype eq 'contig') {
       print_consensus($fh, $uid, $len, $nseq, $lseq, $ftype);
    }

    my $ra_offsets = buildOffsets($lseq);

    # here we parse the individual sequences aligned to the contig
    my($id,$sequence,$asml,$asmr,$rc,$seqleft,$seqright);
    for (my $i = 0; $i <= $#$recs; $i++){
        my ($sid, $sfs, $srecs) = parseCARecord($$recs[$i]);

        if ($sid eq "MPS" and $$sfs{typ} eq "R" ) {
            if ($reverse) {
                ($id,$sequence,$asml,$asmr,$rc,$seqleft,$seqright) =
                   processMPS($sfs,"$len,0", $len);
            } else {
                ($id,$sequence,$asml,$asmr,$rc,$seqleft,$seqright) =
                   processMPS($sfs);
            }
            print_aligned($fh, $seqnames{$id}, $sequence, $asml, $rc,
                          $seqleft, $seqright, $$ra_offsets[$asml],
                          $$ra_offsets[$asmr - 1], $ftype);
        }
        else {
            $tf->logLocal("Record $sid (type $$sfs{typ}) not a unitig" .
                         " sequence, ignoring", 9);
        }
    }
}

######################################################################
# Create the canonical unitig instance name based on class system.
#   input : a number (utg instance), and a base portion of a name.
#   output: a final name for the unitig.
#   so for example, if you have 5 unitigs of the same class,
#   which will begin with "UTG12", they will be assigned the
#   names: UTG12A ($i=0), UTG12B, UTG12C, UTG12D, UTG12E ($i=4).
#   If $i is greater than $ALPHABET_LENGTH, the letters will wrap;
#   ie. UTG12AA (i=26), UTG12AB (i=27) ... UTG12ZY ($i=700) ...
# Written by Dan Kosack.
######################################################################
sub create_name {
    my $i = shift;
    my $name = shift;

    my $ALPHABET_START = 65;                  # ASCII start of alphabet
    my $ALPHABET_LENGTH = 26;

    my @label = ();                           # this goes after $name
    my $num_digits = 1;                       # number of "digits" in label
    my $offset_incr = 0;                      # offset step val., n ^ power
    my $offset_val = 0;                       # offset total
    my $i_copy = $i;                          # copy of i, to modify

    # Calculate the offset and digits.
    my $old_offset = $offset_val;
    while ( $i >= $offset_val ) {
       if ( $offset_incr == 0 ) {
          $offset_incr = 1;
       }
       $old_offset = $offset_val;
       $offset_incr *= $ALPHABET_LENGTH;
       $offset_val += $offset_incr;
       $num_digits++;
    }
    $offset_val = $old_offset;
    $num_digits--;

    # Compute the values for each digit.  Use a dividend / remainder method.
    # Chop off any extraneous lower order material that will only confuse this.
    my $adjusted_i = $i - $offset_val;        # initial value for compute
    my $digit_place;                          # loop controller
    for ( $digit_place = 0; $digit_place < $num_digits; $digit_place++ ) {
       my $i_dividend = int ( $adjusted_i / $ALPHABET_LENGTH );
       my $i_remainder = $adjusted_i - ( $i_dividend * $ALPHABET_LENGTH );
       my $letter_val = chr ($i_remainder + $ALPHABET_START );
       $adjusted_i = $i_dividend;
       unshift @label, $letter_val;           # update with right shifted val
    }
    return $name . join '', @label;
}

######################################################################
# This method outputs a list of degenerate contigs (contigs not in a
# scaffold). It accepts the degenerates list filename and returns
# nothing.
######################################################################
sub outputDegenerates {
    my $degenfname = shift;

    $tf->logLocal("Outputting degenerates information...", 1);
    my $degen_fh = new IO::File ">$degenfname" or
        $tf->bail("Cannot open \"$degenfname\": $!\n");

    foreach my $contig (sort keys %contigs) {
        if ( !exists $scaffold_contigs{$contig} ) {
            print $degen_fh "$contig\n";
        }
    }
    foreach my $contig (sort keys %adjusted_contigs) {
        if ( !exists $scaffold_contigs{$contig} ) {
            print $degen_fh "$contig\n";
        }
    }
    $degen_fh->close();
}

sub calculateUngappedPos($$) {
   my $data = shift;
   my $pos = shift;
   my @data_arr = split(//, $data);
   my $length = scalar @data_arr;
   my $gapped_count = 0;
   my $ungapped_count = 0;
   foreach my $base (@data_arr) {
      if($base ne "-") {
         $ungapped_count++;
      }
      if($gapped_count == $pos) {
         last;
      }
      $gapped_count++;
   }
   return $ungapped_count;
}

######################################################################
# This method outputs the contig/unitig feature report in XML. It
# accepts the report filename. The information is dumped from the
# feats hash.
######################################################################
sub outputFeatures {
    my $featurename = shift;
    my %nu_unitig_count = ();
    $tf->logLocal("Outputting feature information...", 1);
    my $feat_fh = new IO::File ">$featurename" or
        $tf->bail("Cannot open \"$featurename\": $!\n");

    print $feat_fh "$XML_HEADER";
    print $feat_fh "<FeatureSet>\n";

    foreach my $contig_id (sort keys %adjusted_contigs) {
	print $feat_fh "<Contig Id=\"$contig_id\">\n";
       my $ra_utgs = $adjusted_contigs{$contig_id}->{'utgs'};
       ### get consensus sequence for contig
       my($rec,$fields,$recs) = getAsmRecord($contig_id);
       my $contig_data = $$fields{cns};

       for my $utgr (@$ra_utgs) {
          my($utg,$upos) = split /~~/, $utgr;
          my $utype = $unitigs{$utg}->{'type'};
          if ($utype eq 'S') {
            $tf->logLocal("Checking Type S unitig $utg...", 4);
            my $unique = 0;
            my($ustart,$uend) = split /,/, $upos;
            my $type;
            my $name = undef;
            if ($separable_unitigs{$utg} < 2) {  # unique type s unitig
                $tf->logLocal("Checking Unique Type S unitig $utg for " .
                              "gaps...", 4);
                $unique = 1;
            }

            if ($unique) {
               $type = $UNIQUE_FEAT_TYPE;
               $name = $utg;
            }
            else {
               $type = $NONUNIQUE_FEAT_TYPE;
               if (exists $nu_unitig_count{$utg}) {
                    $nu_unitig_count{$utg}++;
               }
               else {
                    $nu_unitig_count{$utg} = 0;
               }
               $name = create_name($nu_unitig_count{$utg},$utg);
            }
            my $unitig_start = $ustart;
            my $unitig_end = $uend;
            my $ungapped_unitig_start = calculateUngappedPos($contig_data, $unitig_start);
            my $ungapped_unitig_end = calculateUngappedPos($contig_data, $unitig_end);
            my $id = "";
            my $class = $utg;
            my $uname = "UTG".$name;
            my $method = $PROGRAM;
            my $assignby = $ENV{'USER'};
            my $comment = "CA_UNITIG ID: $utg";

            print $feat_fh "    <Feature Id=\"$id\" Type=\"$type\" Class=\"$class\" Name=\"$uname\" ".
		           "Method=\"$method\" Assignby=\"$assignby\">\n";
            print $feat_fh "      <Location End5=\"$ungapped_unitig_start\" ".
                            "End3=\"$ungapped_unitig_end\" ".
			    "End5_gapped=\"$unitig_start\" End3_gapped=\"$unitig_end\"/>\n";
            print $feat_fh "      <Comment>$comment</Comment>\n";
            print $feat_fh "    </Feature>\n";
         }
      }
      print $feat_fh "  </Contig>\n";
   }
   print $feat_fh "</FeatureSet>\n";
   $feat_fh->close();
}

######################################################################
# Create the canonical unitig instance name based on class system.
#   input : a number (utg instance), and a base portion of a name.
#   output: a final name for the unitig.
#   so for example, if you have 5 unitigs of the same class,
#   which will begin with "UTG12", they will be assigned the
#   names: UTG12A ($i=0), UTG12B, UTG12C, UTG12D, UTG12E ($i=4).
#   If $i is greater than $ALPHABET_LENGTH, the letters will wrap;
#   ie. UTG12AA (i=26), UTG12AB (i=27) ... UTG12ZY ($i=700) ...
# Written by Dan Kosack.
######################################################################
sub create_name {
    my $i = shift;
    my $name = shift;

    my $ALPHABET_START = 65;                  # ASCII start of alphabet
    my $ALPHABET_LENGTH = 26;

    my @label = ();                           # this goes after $name
    my $num_digits = 1;                       # number of "digits" in label
    my $offset_incr = 0;                      # offset step val., n ^ power
    my $offset_val = 0;                       # offset total
    my $i_copy = $i;                          # copy of i, to modify

    # Calculate the offset and digits.
    my $old_offset = $offset_val;
    while ( $i >= $offset_val ) {
       if ( $offset_incr == 0 ) {
          $offset_incr = 1;
       }
       $old_offset = $offset_val;
       $offset_incr *= $ALPHABET_LENGTH;
       $offset_val += $offset_incr;
       $num_digits++;
    }
    $offset_val = $old_offset;
    $num_digits--;

    # Compute the values for each digit.  Use a dividend / remainder method.
    # Chop off any extraneous lower order material that will only confuse this.
    my $adjusted_i = $i - $offset_val;        # initial value for compute
    my $digit_place;                          # loop controller
    for ( $digit_place = 0; $digit_place < $num_digits; $digit_place++ ) {
       my $i_dividend = int ( $adjusted_i / $ALPHABET_LENGTH );
       my $i_remainder = $adjusted_i - ( $i_dividend * $ALPHABET_LENGTH );
       my $letter_val = chr ($i_remainder + $ALPHABET_START );
       $adjusted_i = $i_dividend;
       unshift @label, $letter_val;           # update with right shifted val
    }
    return $name . join '', @label;
}

######################################################################
# This method generates the unitig report. It accepts the report
# filename. The information is dumped from the contigs and unitigs
# hashes.
######################################################################
sub generateUnitigReport($$) {
    my $unitigname = shift;
    my $unitig_config = shift;
    my %unitig_hash = ();
    my $unitig_fh = new IO::File ">$unitigname" or
        $tf->bail("Cannot open \"$unitigname\": $!\n");

    my $config_fh = new IO::File ">$unitig_config" or
       $tf->bail("Cannot open \"$unitig_config\": $!\n");

    $tf->logLocal("Generating report of contig/unitig relationships.", 1);
    print $unitig_fh '*' x 70 . "\n";
    print $unitig_fh "Unitig Report\n";
    print $unitig_fh '*' x 70 . "\n";

        print $unitig_fh qq~
-----------------------------
Contigs
-----------------------------
~;
    foreach my $k (sort keys %contigs) {
        my $v = $contigs{$k};
        print $unitig_fh qq~Contig: $k
Length: $v->{'len'}
Sequences: $v->{'nseq'}
Unitig(s):
~;

    foreach my $utgr (@{$v->{'utgs'}}) {
       my($utg,$upos) = split /~~/, $utgr;
       print $unitig_fh "           $utg (Typ: $unitigs{$utg}->{'type'}".
                        " Seqs: $unitigs{$utg}->{'nseq'} Len: " .
                        "$unitigs{$utg}->{'len'} Pos: $upos)\n";
       my $contig_len = $v->{'len'};
       my $unitig_count = -1;
       if(defined $unitig_hash{$utg}) {
	  $unitig_count = $unitig_hash{$utg};
          $unitig_count++;
          $unitig_hash{$utg} = $unitig_count;
       }
       else {
          $unitig_count++;
          $unitig_hash{$utg} = $unitig_count;
       }
       printConfigInfo($k, $contig_len, $utg, $upos, $config_fh,
                       $unitig_count);
    }
print $unitig_fh qq~
-----------------------------
~;
    }

print $unitig_fh qq~
-----------------------------
Adjusted Contigs
-----------------------------
~;
    foreach my $k (sort keys %adjusted_contigs) {
        my $v = $adjusted_contigs{$k};
        print $unitig_fh qq~Contig: $k
Sequences: $v->{'nseq'}
Unitig(s):
~;
        my @separables = ();
        foreach my $utgr (@{$v->{'utgs'}}) {
            my($utg,$upos) = split /~~/, $utgr;
            my $utype = $unitigs{$utg}->{'type'};
            print $unitig_fh "           $utg (Typ: $utype".
                     " Seqs: $unitigs{$utg}->{'nseq'} Len: " .
                     "$unitigs{$utg}->{'len'} Pos: $upos)\n";
            if ($utype eq 'S') {
                push(@separables, $utg);
            }
        }
        print $unitig_fh "Separable (Type S) unitigs:\n";

        foreach my $sep_utg (@separables) {
            print $unitig_fh "                            " .
                             "Unitig $sep_utg ";
            if ($separable_unitigs{$sep_utg} > 1) {
                print $unitig_fh "[Non-unique]\n";
            } else {
                print $unitig_fh "[Unique]\n";
            }
        }
        print $unitig_fh qq~
-----------------------------
~;
    }

    print $unitig_fh qq~
-----------------------------
Unitigs Promoted to Contigs
-----------------------------
~;
    foreach my $k (sort keys %separable_unitigs) {
       if($separable_unitigs{$k} == 1) {
          my $in_contig = '';
          foreach my $ctgr (@{$unitigs{$k}->{'in_contig'}}) {
             my $ctg = (split /~~/, $ctgr)[0];
             $in_contig .= "$ctg ";
          }

          print $unitig_fh qq~Contig: $k
Length: $unitigs{$k}->{'len'}
Seqs: $unitigs{$k}->{'nseq'}
In contig(s): $in_contig

~;
       }
    }
    print $unitig_fh qq~
-----------------------------
Totals
-----------------------------
~;
    my $total_not_in_contig = 0;
    my @not_in_contig = ();
    my $total_sep_unitigs = scalar keys %separable_unitigs;
    while (my($k,$v) = each %unitigs) {
        unless (scalar @{$v->{'in_contig'}} > 0) {
            $total_not_in_contig++;
            push(@not_in_contig, $k);
        }
    }
    my $total_contigs = keys %contigs;
    my $total_adj_contigs = keys %adjusted_contigs;
    my $total_unitigs = keys %unitigs;
    print $unitig_fh qq~
Total Contigs: $total_contigs
Total Contigs with Separable Unitigs: $total_adj_contigs
Total Unitigs: $total_unitigs
Total Separable Unitigs: $total_sep_unitigs
Total Unitigs not in a contig: $total_not_in_contig
~;
    if ($total_not_in_contig) {
        print $unitig_fh "Unitigs Not in a contig: @not_in_contig\n";
    }
$unitig_fh->close();
}

sub printConfigInfo {
    my $contig = shift;
    my $contig_len = shift;
    my $utg = shift;
    my $upos = shift;
    my $unitig_fh = shift;
    my $unitigcount = shift;

    my $ulen = $unitigs{$utg}->{'len'};
    my $u_gaps = $unitigs{$utg}->{'gaps'};
    my $contig_type = undef;
    if (exists $scaffold_contigs{$contig}) {
       $contig_type = "CA_CONTIG";
    }
    else {
       $contig_type = "CA_DEGEN";
    }

    # printing the unitig information in the unitig config file
    print $unitig_fh "[unitig\_$unitigcount\_$utg]\n";
    print $unitig_fh "contig_id=$contig\n";
    print $unitig_fh "contig_length=$contig_len\n";
    print $unitig_fh "contig_type=$contig_type\n";

    my @aligns = split(/,/,$upos);
    my $start = $aligns[0];
    my $end =  $aligns[1];

    print $unitig_fh "contig_alignstart=$start\n";
    print $unitig_fh "contig_alignend=$end\n";

    $u_gaps =~ s/\n/ /g;
    my @gaps = split(/ /,$u_gaps);
    my $gaps = join ",", @gaps;

    print $unitig_fh "unitig_length=$ulen\n";
    print $unitig_fh "unitig_gaps=$gaps\n";

    if(defined $separable_unitigs{$utg}) {

       if ($separable_unitigs{$utg} > 1) {
          my $unitig_cons = getUnitigSequence($utg, $start, $end);
          print $unitig_fh "unitig_cons=$unitig_cons\n";
          print $unitig_fh "unitig_type=NON_UNIQUE_SURROGATE\n";
          print $unitig_fh "unitig_promotion=not promoted\n";
       }
       else {
          my $unitig_cons = getUnitigSequence($utg, $start, $end);
          print $unitig_fh "unitig_cons=$unitig_cons\n";
          print $unitig_fh "unitig_type=UNIQUE_SURROGATE\n";
          print $unitig_fh "unitig_promotion=promoted\n";
       }
    }
    else {
       print $unitig_fh "unitig_type=NON_SURROGATE\n";
       print $unitig_fh "unitig_promotion=promoted\n";
    }
    print $unitig_fh "\n";
}

sub getUnitigSequence {
   my $unitig_id = shift;
   my $start = shift;
   my $end = shift;
   my($rec,$fields,$recs) = getAsmRecord($unitig_id);
   my $unitig = $$fields{cns};
   $unitig =~ s/\n//g;
   $unitig =~ s/\s*//g;
   if ($start > $end) {
      my $rev_unitig = reverseComplement($unitig);
      $unitig = $rev_unitig;
   }
   return $unitig;
}

######################################################################
# This method swaps two numbers. It accepts a reference to a scalar
# containing the first number and a reference to a scalar containing
# the second number.
######################################################################
sub swap {
    my $r_a = shift;
    my $r_b = shift;
    my $tmp;

    $tmp = $$r_a;
    $$r_a = $$r_b;
    $$r_b = $tmp;
}

######################################################################
# This method calculates the redundancy information for a contig.
# It accepts a contig filename.
######################################################################
sub calculateContigInfo {
    my $contigname = shift;
    $tf->logLocal("Calculating redundancy info for $contigname", 4);
    my @lines = `grep \"#\" $contigname`;

    my $contig_line = undef;
    my %contig_info = ();
    my %contig_ints = ();
    my $prev_asm_right = undef;
    my $contig_id = undef;
    my $contig_len = undef;
    my $total_seq_len = 0;
    my $asm_left = undef;
    my $asm_right = undef;
    my $num_bases = undef;
    my $contig_num_bases = undef;
    my $max_asmr = 0;

    foreach $contig_line (@lines) {

       if($contig_line =~ /$CONTIG_HEADER_PARSE/) {
          if(defined $contig_id) {
             if($max_asmr < $contig_num_bases) {
                my $non_cov_len = $contig_num_bases - $max_asmr;
                $contig_len -= $non_cov_len;
             }
             my $red = 0;
             if($contig_len != 0) {
                $red = ($total_seq_len/$contig_len);
	     }
             my $redundancy = sprintf("%.2f",$red);
             $total_seq_len = 0;
             $max_asmr = 0;
             $prev_asm_right = undef;
             $contig_num_bases = undef;
             $tasm_subs{$contig_id} = $redundancy;
          }
          $contig_id = $1;
          $contig_len = $3;
          $contig_num_bases = $contig_len;
       }
       elsif($contig_line =~ /$CONTIG_SEQ_HEADER_PARSE/) {
          my $seq_id = $1;
          my $offset = $2;
          $num_bases = $4;
          $asm_left = $offset + 1;
          $asm_right = $offset + $num_bases;
          if($asm_right > $max_asmr) {
             $max_asmr = $asm_right;
          }

          $total_seq_len += (($asm_right - $asm_left) + 1);

          if(defined $prev_asm_right) {
             if(($asm_left > $prev_asm_right) &&
                ($prev_asm_right > $max_asmr)) {
	        my $non_cov_len = ($asm_left - $prev_asm_right) - 1;
                $contig_len -= $non_cov_len;
	     }
          }
          else { # first sequence
	     if($asm_left > 1) {
                my $non_cov_len = ($asm_left - 1);

                $contig_len -= $non_cov_len;
             }
          }
          $prev_asm_right = $asm_right;
      }
   }

   if(defined $contig_id) {

      if($max_asmr < $contig_num_bases) {
         my $non_cov_len = $contig_num_bases - $max_asmr;
         $contig_len -= $non_cov_len;
      }
      my $red = 0;
      if($contig_len != 0) {
         $red = ($total_seq_len/$contig_len);
      }
      my $redundancy = sprintf("%.2f",$red);
      $total_seq_len = 0;
      $prev_asm_right = undef;
      $contig_num_bases = undef;
      $tasm_subs{$contig_id} = $redundancy;
   }
}

######################################################################
# This method processes the MPS records found in CCO and UTG records.
# It accepts a reference to a hash containing the fields of the
# record, an optional global unitig position, and an optional string
# of contig gaps to add. It returns the sequence name, the sequence,
# the asm left and right, the reverse complement flag, and the seq
# left and right.
######################################################################
sub processMPS {
    my $sfs = shift;
    my $upos = shift || undef;
    my $ulen = shift || undef;
    my $cgaps = shift || undef;
    my $id = getCAId($$sfs{mid});
    my $asms = $$sfs{pos};

    $asms =~ /(\d+),(\d+)/;
    if (! defined $1){
        $tf->bail("Badly formed position record: $$sfs{pos}");
    }
    my $asml = $1;
    my $asmr = $2;
    $tf->logLocal("Got asm_left = $asml, asm_right = $asmr for $id " .
                  "($seqnames{$id})...", 9);
    my $gapno = $$sfs{dln};
    $tf->logLocal("Got $gapno gaps for $id ($seqnames{$id})...", 9);
    my $gaps = $$sfs{del};
    my @gaps = split(' ', $gaps);
    if ($gapno != $#gaps + 1){
        $tf->logError("Warning: $gapno != " . $#gaps + 1);
    }
    my $sequence = get_seq($frag_fh, $id);
    my @lines = split('\n', $sequence);
    $sequence = join('', @lines);
    $sequence = uc($sequence);
    $sequence = substr($sequence, $seql{$id}, $seqr{$id} - $seql{$id});

    my $seqleft = $seql{$id} + 1;
    my $seqright = $seqr{$id};
    my $rc = "";

    # sequence is reverse complemented
    if ($asml > $asmr) {
        $sequence  = reverseComplement($sequence);
        my $tmp = $asmr;
        $asmr = $asml;
        $asml = $tmp;
        $tmp = $seqright;
        $seqright = $seqleft;
        $seqleft = $tmp;
        $rc = "RC";
    }

    # add in gaps in sequence
    my $gaps_added = 0;
    ($sequence, $gaps_added) = addUGaps($sequence, \@gaps);
    $tf->logLocal("Added $gaps_added gaps for $id ($seqnames{$id})...", 9);

    # adjust asmr in case we have gaps at the very end of the sequence
    # that we've ignored
    if ($asmr - $asml > length($sequence)) {
        $asmr = $asml + length($sequence) - 1;
        $tf->logLocal("adjusting $asmr to $asmr", 1);
    }

    if ($asmr - $asml < length($sequence)) {
        $tf->logError("Looks like a major bug: Sequence $id, " .
                      ($asmr - $asml) . " != " . length($sequence) . " !!!");
    }

    # if a contig position is passed,
    # aligning unitig seqs to the consensus:
    #   shifting asml and asmr
    #   reverse complementing if necessary
    #   adding in gaps from unitig consensus
    if (defined $upos) {
        my($global_start,$global_end);
        my $align = '';
        my($ul,$ur) = split /,/, $upos;
        $tf->logLocal("Aligning unitig seq $id ($seqnames{$id}) to " .
                      "contig, global position $ul to $ur", 4);
        if ($ul < $ur) {
            $global_start = $ul;
            $global_end = $ur;
        } else {
            $global_start = $ur;
            $global_end = $ul;
            $align = 'reverse';
        }

        $tf->logLocal("Original asm_left = $asml, asm_right = $asmr ".
                      "for $id ($seqnames{$id})...", 4);
        # adjust asml, asmr if unitig is reversed
        my $local_offset = $asml < $asmr ? $asml : $asmr;
        if ($align eq 'reverse') {
            $sequence  = reverseComplement($sequence);
            $rc = $rc eq "" ? "RC" : "";
            $asml = $ulen - $local_offset - length($sequence);
            $asmr = $asml + length($sequence);
            my $tmp = $seqright;
            $seqright = $seqleft;
            $seqleft = $tmp;
            $tf->logLocal("Unitig reversed: local asm_left = $asml, " .
                          "local asm_right = $asmr", 4);
        }

        # add in any extra gaps needed to align the unitig consensus with
        # the contig
        if (defined $cgaps) {
            my @cgaps = split(' ', $cgaps);
            my @sgaps = ();
            my $prev_gaps = 0;
            foreach my $cgap (@cgaps) {
                if ($cgap < $asml) {
                    $prev_gaps++;
                } elsif ($cgap >= $asml and $cgap <= $asmr) {
                    my $sgap = $cgap - $asml;
                    $tf->logLocal("Adding extra gap at $sgap for seq " .
                                  "$id ($seqnames{$id}) to align unitig " .
                                  "with the contig", 4);
                    push(@sgaps, $sgap);
                }
            }
            # adjust local asml, asmr for gaps added before seq start
            $asml += $prev_gaps;
            $asmr += $prev_gaps;
            $tf->logLocal("For seq $id ($seqnames{$id}), $prev_gaps gaps " .
                          "were inserted prior to seq start; " .
                          "local asml = $asml, local asmr = $asmr", 4);

            # add in gaps
            my $num_gaps = scalar @sgaps;
            if ($num_gaps > 0) {
                my $gaps_added = 0;
                ($sequence, $gaps_added) = addGaps($sequence, \@sgaps);
                $asmr += $gaps_added;
                $tf->logLocal("$gaps_added extra gaps were added to seq $id ".
                              "($seqnames{$id}), local asmr = $asmr", 4);
            }
        }

        # adjusting asml and asmr to global positions
        $asml += $global_start;
        $asmr += $global_start;
        $tf->logLocal("Adjusted asm_left = $asml, asm_right = $asmr " .
                      "for $id ($seqnames{$id})...", 4);
    }
    return ($id, $sequence, $asml, $asmr, $rc, $seqleft, $seqright);
}

######################################################################
# This method adds gaps into a sequence (using gapped coordinates).
# It accepts a sequence and a reference to an array containing the
# gaps to add. It returns the sequence with the gaps added.
# Gaps are counted in the following manner:
#
# A A A T G C - A G C
# 1 2 3 4 5 6 6 7 8 9
######################################################################
sub addGaps {
    my $sequence = shift;
    my $ra_gaps = shift;
    my @gseq = ();
    my $added = 0;
    if (scalar @$ra_gaps > 0) {
        @gseq = split //, $sequence;
        my $ungapped_len = scalar @gseq;
        foreach my $gap (@$ra_gaps) {
            if ($gap >= $ungapped_len) { last; }
            splice(@gseq, $gap, 0, '-');
            $ungapped_len++;
            $added++;
        }
        $sequence = join '', @gseq;
    }
    return $sequence, $added;
}

######################################################################
# This method adds gaps into a sequence (using ungapped coordinates).
# It accepts a sequence and a reference to an array containing the
# gaps to add. It returns the sequence with the gaps added, the number
# of gaps inserted, and the gapped positions of the inserted gaps.
# Gaps are counted in the following manner:
#
# A A A T G C - A G C
# 1 2 3 4 5 6 6 7 8 9
######################################################################
sub addUGaps {
    my $sequence = shift;
    my $ra_gaps = shift;
    my @gseq = ();
    my @gap_pos = ();
    my($shift,$ungapped_len);
    if (scalar @$ra_gaps > 0) {
        my $ra_curr_gaps = findUGaps($sequence);
        @gseq = split //, $sequence;
        $ungapped_len = scalar @gseq;
        $shift = 0;
        foreach my $gap (@$ra_gaps) {
            my $ra_prev_gaps = getPrevGaps($ra_curr_gaps, $gap);
            my $num_prev_gaps = scalar @$ra_prev_gaps;
            my $insertpos = $gap + $num_prev_gaps + $shift;
            if ($insertpos >= $ungapped_len) { last; }
            push(@gap_pos, $insertpos);
            splice(@gseq, $insertpos, 0, '-');
            $shift++;
            $ungapped_len++;
        }
        $sequence = join '', @gseq;
    }
    return $sequence, $shift, \@gap_pos;
}

######################################################################
# This method finds the number of gaps in a sequence before the
# current position. It accepts a reference to an array of current gaps
# and a position. It returns a reference to an array containing the
# the previous gaps.
######################################################################
sub getPrevGaps {
    my $ra_gaps = shift;
    my $index = shift;
    my @prev_gaps = ();

    foreach my $gap (@$ra_gaps) {
        if ($gap >= $index) { last; }
        push(@prev_gaps, $gap);
    }
    return \@prev_gaps;
}

######################################################################
# This method finds the (ungapped) positions of gaps in a sequence.
# It accepts a sequence. It returns a reference to an array containing
# the gaps positions.
######################################################################
sub findUGaps {
    my $seq = shift;

    my $regex = "-";
    my $shift = 0;
    my @index = ();
    while ($seq =~ m/$regex/g) {
        my $offset = length($`) - $shift;
        push(@index, $offset);
        $shift++;
    }
    return \@index;
}

######################################################################
# This method prints a contig consensus in either .contig or .asm
# format to the .contig or .tasm file. It accepts a filehandle, a
# sequence name, the length of the sequence, the number of sequences,
# the sequence, the type (contig|asm) which describes the format, an
# optional comment, and a flag (0|1) which signifies if this is the
# first record printed.
######################################################################
sub print_consensus {
    my $file = shift;
    my $id = shift;
    my $len = shift;
    my $nseq = shift;
    my $sequence = shift;
    my $how = shift;
    my $type = shift or 0;
    my $first = shift or 0;

    if ($how eq "contig") {
        print $file "\#\#$id $nseq $len bases, 00000000 checksum.\n";
        print_sequence($file, $sequence);
    } elsif ($how eq "asm") {
        my $strip = $sequence;
        $strip =~ s/-//g;

        my $quality = "";
        my $redundancy = "";
        if(defined $tasm_subs{$id}) {
	   $redundancy = $tasm_subs{$id};
        }
        $quality = "0x";
        for (my $i = 0; $i < length($sequence); $i++) {
           $quality .= "06";
        }

        if ( ! $first ) {
            print $file "|\n";
        }
        print $file qq~sequence\t$strip
lsequence\t$sequence
asmbl_id\t$id
seq_id\t
com_name\t
type\t
method\tCelera Assembler
ed_status\t
redundancy\t$redundancy
perc_N\t
seq#\t$nseq
full_cds\t
cds_start\t
cds_end\t
ed_pn\t$who
ed_date\t$date
comment\t$type ID: $id
frameshift\t
ca_contig_id\t$id
mod_date\t$date
~;
        return;
    } elsif ($how eq "fasta") {
        $sequence =~ s/-//g;  # get rid of all gaps
        my $fasta_len = length($sequence);
        print $file ">$id $nseq $fasta_len bases\n";
        print_sequence($file, $sequence);
    }
}

######################################################################
# This method prints an aligned sequence to a contig or a unitig to
# the .contig or .tasm file. It accepts a filehandle, the sequence
# name, the offset, the reverse complement flag, the seq left and
# right, the asm left and right, and the type (contig|asm) which
# describes the print format.
######################################################################
sub print_aligned {
    my($file, $name, $seq, $offset, $rc,
       $seqleft, $seqright, $asml, $asmr, $type) = @_;

    if ($type eq "contig"){
        print $file "\#$name($offset) [$rc] ",
          length($seq),
          " bases, 00000000 checksum. {$seqleft $seqright} <$asml $asmr>\n";
        print_sequence($file, $seq);
    } elsif ($type eq "asm"){
        print $file qq~
seq_name\t$name
asm_lend\t$asml
asm_rend\t$asmr
seq_lend\t$seqleft
seq_rend\t$seqright
best\t
comment\t
db\t
offset\t$offset
lsequence\t$seq
~;
    }
}

######################################################################
# This method prints a sequence in 60 character per line format. It
# accepts a filehandle and a sequence.
######################################################################
sub print_sequence {
    my $file = shift;
    my $seqs = shift;

    $seqs =~ s/((.){60})/$1\n/g;
    $seqs =~ s/\n$//;
    print $file "$seqs\n";
}

######################################################################
# This method gets a sequence fragment from the frg file. It accepts
# a filehandle and an id. The id is a key into a hash containing the
# position in the frg file. It returns the sequence.
######################################################################
sub get_seq {
    my $file = shift;
    my $id = shift;
    my $pos = $seqpos{$id};
    seek $file, $seqpos{$id}, 0; # seek set
    my $record = getCARecord($file);
    if (! defined $record){
        $tf->bail("Error: trouble getting seq record from $fragfname");
    }

    my ($rec, $fields, $recs) = parseCARecord($record);
    if ($rec ne "FRG"){
        $tf->bail("Error in get_seq, expecting frg record, got $rec.");
    }
    if ($$fields{acc} != $id){
        $tf->bail("Error in get_seq, expecting $id, got $$fields{acc}");
    }
    return $$fields{seq};
}

######################################################################
# This method seeks to a certain position in the file to be read,
# reads a buffer of a specified length, and prints the buffer out to
# the file to be written. It returns nothing.
######################################################################
sub printBufferToFile {
    my $read_fh = shift;
    my $start = shift;
    my $end = shift;
    my $write_fh = shift;
    my $buffer = "";
    my $length = $end - $start + 1;

    $read_fh->seek($start - 1, 0);
    $read_fh->read($buffer, $length);

    print $write_fh $buffer;
}
