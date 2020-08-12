
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

package canu::Report;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(generateReport addToReport getFromReport);

use strict;
use warnings "all";
no  warnings "uninitialized";

#use File::Copy;
#use File::Path 2.08 qw(make_path remove_tree);

#use canu::Defaults;
use canu::Grid_Cloud;


#  Holds the report we have so far (loaded from disk) and is then populated with
#  any new results and written back to disk.
#
#  For a rather technical reason, some keys will exist in the hash (*Meryl, for example)
#  but have an undef value.
#
my %report;



#  Parse an existing report to load the data present in it.

sub loadReport ($) {
    my $asm = shift @_;
    my $tag;
    my $rpt;

    fetchFile("$asm.report");

    return  if (! -e "$asm.report");

    #  Parse the file, loading each section into the proper key.

    open(F, "< $asm.report") or caExit("can't open '$asm.report' for reading: $!", undef);
    while (<F>) {
        if      (m/^\[CORRECTION/) {
            $tag = "cor";
        } elsif (m/^\[TRIMMING/) {
            $tag = "obt";
        } elsif (m/^\[UNITIGGING/) {
            $tag = "utg";
        }

        if      (m/READS\]$/)        {  $rpt = "${tag}SeqStore";  $report{$rpt} = undef;
        } elsif (m/MERS\]$/)         {  $rpt = "${tag}Meryl";     $report{$rpt} = undef;

        } elsif (m/FILTERING\]$/)    {  $rpt = "corFilter";       $report{$rpt} = undef;
        } elsif (m/LAYOUT\]$/)       {  $rpt = "corLayout";          $report{$rpt} = undef;
        } elsif (m/CORRECTIONS\]$/)  {  $rpt = "corrections";     $report{$rpt} = undef;

        } elsif (m/TRIMMING\]$/)     {  $rpt = "trimming";        $report{$rpt} = undef;
        } elsif (m/SPLITTING\]$/)    {  $rpt = "splitting";       $report{$rpt} = undef;

        } elsif (m/OVERLAPS\]$/)     {  $rpt = "overlaps";        $report{$rpt} = undef;
        } elsif (m/ADJUSTMENT\]$/)   {  $rpt = "adjustments";     $report{$rpt} = undef;
        } elsif (m/ERROR RATES\]$/)  {  $rpt = "error rates";     $report{$rpt} = undef;
        } elsif (m/UNITIGS\]$/)      {  $rpt = "unitigs";         $report{$rpt} = undef;
        } elsif (m/CONTIGS\]$/)      {  $rpt = "contigs";         $report{$rpt} = undef;
        } elsif (m/CONSENSUS\]$/)    {  $rpt = "consensus";       $report{$rpt} = undef;

        } else {
            $report{$rpt} .= $_;
        }
    }
    close(F);

    #  Ignore blank lines at the start and end.

    foreach my $k (keys %report) {
        $report{$k} =~ s/^\s+//;
        $report{$k} =~ s/\s+$//;
        $report{$k} .= "\n";
    }
}



sub saveReportItem ($$) {
    my $title = shift @_;
    my $item  = shift @_;

    print F "[$title]\n$item\n"   if (defined($item));
}



sub saveReport ($) {
    my $asm = shift @_;
    my $tag;

    open(F, "> $asm.report.new") or caExit("can't open '$asm.report.new' for writing: $!", undef);

    saveReportItem("CORRECTION/READS",       $report{"corSeqStore"});
    saveReportItem("CORRECTION/MERS",        $report{"corMeryl"});
    saveReportItem("CORRECTION/FILTERING",   $report{"corFilter"});
    saveReportItem("CORRECTION/LAYOUT",      $report{"corLayout"});
    saveReportItem("CORRECTION/CORRECTIONS", $report{"corrections"});

    saveReportItem("TRIMMING/READS",         $report{"obtSeqStore"});
    saveReportItem("TRIMMING/MERS",          $report{"obtMeryl"});
    saveReportItem("TRIMMING/TRIMMING",      $report{"trimming"});
    saveReportItem("TRIMMING/SPLITTING",     $report{"splitting"});

    saveReportItem("UNITIGGING/READS",       $report{"utgSeqStore"});
    saveReportItem("UNITIGGING/MERS",        $report{"utgMeryl"});
    saveReportItem("UNITIGGING/OVERLAPS",    $report{"overlaps"});
    saveReportItem("UNITIGGING/ADJUSTMENT",  $report{"adjustments"});
    saveReportItem("UNITIGGING/ERROR RATES", $report{"error rates"});
    saveReportItem("UNITIGGING/UNITIGS",     $report{"unitigs"});
    saveReportItem("UNITIGGING/CONTIGS",     $report{"contigs"});
    saveReportItem("UNITIGGING/CONSENSUS",   $report{"consensus"});

    close(F);

    my $diff = -1;

    if (-e "$asm.report") {
        open(OLD, "< $asm.report");
        open(NEW, "< $asm.report.new");

        $diff = 0;

        while (!eof(OLD) || !eof(NEW)) {
            my $old = <OLD>;
            my $new = <NEW>;

            if ($old ne $new) {
                $diff++;
            }
        }

        close(OLD);
        close(NEW);
    }


    if ($diff < 0) {
        #print STDERR "-- New report created.\n";

        rename "$asm.report.new", "$asm.report";
        stashFile("$asm.report");
    }

    elsif ($diff == 0) {
        #print STDERR "-- No change in report.\n";

        unlink "$asm.report.new";
    }

    else {
        #print STDERR "-- Report changed.\n";

        unlink "$asm.report";
        rename "$asm.report.new", "$asm.report";

        stashFile("$asm.report");
    }
}



sub generateReport ($) {
    my $asm = shift @_;

    loadReport($asm);
    saveReport($asm);
}



sub addToReport ($$) {
    my $item = shift @_;
    my $text = shift @_;

    print STDERR $text;     #  Client code is GREATLY simplified if this dumps to the screen too.

    $report{$item} = $text;
}



sub getFromReport ($) {
    return($report{$_[0]});
}


1;
