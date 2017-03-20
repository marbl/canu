
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-NOV-08
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Report;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(generateReport addToReport getFromReport);

use strict;

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

        if      (m/READS\]$/)        {  $rpt = "${tag}GkpStore";  $report{$rpt} = undef;
        } elsif (m/MERS\]$/)         {  $rpt = "${tag}Meryl";     $report{$rpt} = undef;

        } elsif (m/FILTERING\]$/)    {  $rpt = "filtering";       $report{$rpt} = undef;
        } elsif (m/CORRECTIONS\]$/)  {  $rpt = "corrections";     $report{$rpt} = undef;

        } elsif (m/TRIMMING\]$/)     {  $rpt = "trimming";        $report{$rpt} = undef;
        } elsif (m/SPLITTING\]$/)    {  $rpt = "splitting";       $report{$rpt} = undef;

        } elsif (m/OVERLAPS\]$/)     {  $rpt = "overlaps";        $report{$rpt} = undef;
        } elsif (m/ADJUSTMENT\]$/)   {  $rpt = "adjustments";     $report{$rpt} = undef;
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

    open(F, "> $asm.report") or caExit("can't open '$asm.report' for writing: $!", undef);
    
    saveReportItem("CORRECTION/READS",       $report{"corGkpStore"});
    saveReportItem("CORRECTION/MERS",        $report{"corMeryl"});
    saveReportItem("CORRECTION/FILTERING",   $report{"correctionFiltering"});
    saveReportItem("CORRECTION/CORRECTIONS", $report{"corrections"});

    saveReportItem("TRIMMING/READS",         $report{"obtGkpStore"});
    saveReportItem("TRIMMING/MERS",          $report{"obtMeryl"});
    saveReportItem("TRIMMING/TRIMMING",      $report{"trimming"});
    saveReportItem("TRIMMING/SPLITTING",     $report{"splitting"});

    saveReportItem("UNITIGGING/READS",       $report{"utgGkpStore"});
    saveReportItem("UNITIGGING/MERS",        $report{"utgMeryl"});
    saveReportItem("UNITIGGING/OVERLAPS",    $report{"overlaps"});
    saveReportItem("UNITIGGING/ADJUSTMENT",  $report{"adjustments"});
    saveReportItem("UNITIGGING/UNITIGS",     $report{"unitigs"});
    saveReportItem("UNITIGGING/CONTIGS",     $report{"contigs"});
    saveReportItem("UNITIGGING/CONSENSUS",   $report{"consensus"});

    close(F);

    stashFile("$asm.report");
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
