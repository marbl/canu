
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
 #  This file is derived from:
 #
 #    src/pipelines/ca3g/Gatekeeper.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2016-NOV-06
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::HTML;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(buildHTML);

use strict;

use canu::Defaults;
use canu::Execution;


sub buildGatekeeperHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Input Reads</h2>\n";
    push @$body, "\n";

    if (! -e "$wrk/$asm.gkpStore/load.dat") {
        push @$body, "<p>None loaded.</p>\n";
        return;
    }

    push @$body, "<table>\n";

    open(F, "< $wrk/$asm.gkpStore/load.dat") or caExit("can't open '$wrk/$asm.gkpStore/load.dat' for reading: $!", undef);
    while (<F>) {

        #  nam blocks show up once per file.
        if (m/^nam\s(\d+)\s(.*)$/) {
            my $idx  = $1;
            my $file = $2;

            push @$body, "<tr id='gkpload$idx'><td colspan='2'>$file</td></tr>\n";

            push @$scripts, "document.getElementById('gkpload$idx').onclick = toggleTable;\n";
            push @$scripts, "document.getElementById('gkpload$idx').style   = 'cursor: pointer;';\n";
        }

        #  lib blocks show up once per file, all paramters are on the same line
        elsif (m/^lib\s/) {
            my @libs = split '\s+', $_;
            my ($param, $np, $var, $val);

            $param = shift @libs;  #  Throw out the first 'lib' word.
            $np    = scalar(@libs);
            $param = shift @libs;  #  First thing we want to report.

            #  First row needs to have a spanning cell for the 'parameters'.
            ($var, $val) = split '=', $param;
            push @$body, "<tr class='details'><td rowspan='$np'>Parameters</td><td>$var = $val</td></tr>\n";

            #  Remaining rows just have var=val.
            foreach $param (@libs) {
                ($var, $val) = split '=', $param;
                push @$body, "<tr class='details'><td>$var = $val</td></tr>\n";
            }
        }

        #  dat blocks show up once per file, and are the last block emitted for a file
        elsif (m/^dat\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my $nLOADEDA  = $1;
            my $bLOADEDA  = $2;
            my $nSKIPPEDA = $3;
            my $bSKIPPEDA = $4;
            my $nLOADEDQ  = $5;
            my $bLOADEDQ  = $6;
            my $nSKIPPEDQ = $7;
            my $bSKIPPEDQ = $8;
            my $nWARNS    = $9;

            push @$body, "<tr class='details'><td rowspan='2'>FASTA</td><td>$nLOADEDA reads ($bLOADEDA bp)</td></tr>\n",;
            push @$body, "<tr class='details'><td>$nSKIPPEDA reads ($bSKIPPEDA bp) were short and not loaded</td></tr>\n";

            push @$body, "<tr class='details'><td rowspan='2'>FASTQ</td><td>$nLOADEDQ reads ($bLOADEDQ bp)</td></tr>\n";
            push @$body, "<tr class='details'><td>$nSKIPPEDQ reads ($bSKIPPEDQ bp) were short and not loaded</td></tr>\n";

            my $nl = $nLOADEDA  + $nLOADEDQ;
            my $bl = $bLOADEDA  + $bLOADEDQ;
            my $ns = $nSKIPPEDA + $nSKIPPEDQ;
            my $bs = $bSKIPPEDA + $bSKIPPEDQ;

            push @$body, "<tr><td colspan='2'>$nl reads ($bl bp) loaded, $ns reads ($bs bp) skipped, $nWARNS warnings</td></tr>\n";
        }
        
        #  the sum block shows up excatly once, a summary of all the reads loaded
        elsif (m/^sum\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)$/) {
            my $nLOADED  = $1;
            my $bLOADED  = $2;
            my $nSKIPPED = $3;
            my $bSKIPPED = $4;
            my $nWARNS   = $5;

            push @$body, "</table>\n";
            push @$body, "\n";
            push @$body, "<h2>Final Store</h2>\n";
            push @$body, "\n";
            push @$body, "<table>\n";
            push @$body, "<tr><td colspan='2'>$wrk/$asm.gkpStore</td></tr>\n";
            push @$body, "<tr><td>readsLoaded</td><td>$nLOADED reads ($bLOADED bp)</td></tr>\n";
            push @$body, "<tr><td>readsSkipped</td><td>$nSKIPPED reads ($bSKIPPED bp) (read was too short)</td></tr>\n";
            push @$body, "<tr><td>warnings</td><td>$nWARNS warnings (invalid base or quality value)</td></tr>\n";
            push @$body, "</table>\n";
 
        } else {
            caExit("failed to read '$wrk/$asm.gkpStore/load.log': invalid format", undef);
        }
    }
    close(F);

    push @$body, "<h3>Read Length Histogram</h3>\n";
    push @$body, "<a href='$wrk/$asm.gkpStore/readlengths.lg.png'><img src='$wrk/$asm.gkpStore/readlengths.sm.png'></img></a>\n";
}


sub buildMerylHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>k-Mer Counts</h2>\n";
    push @$body, "\n";

    if (! -d "$wrk/0-mercounts") {
        push @$body, "<p>Stage not computed. ($wrk/0-mercounts)</p>\n";
        return;
    }

    my %merSizes;

    open(F, "ls $wrk/0-mercounts/ |") or caExit("can't find files in '$wrk/0-mercounts': $!", undef);
    while (<F>) {
        if (m/\.ms(\d+)\./) {
            $merSizes{$1}++;
        }
    }
    close(F);

    foreach my $ms (keys %merSizes) {
        my $numTotal    = 0;
        my $numDistinct = 0;
        my $numUnique   = 0;
        my $largest     = 0;

        if (-e "$wrk/0-mercounts/$asm.ms$ms.histogram.info") {
            open(F, "<  $wrk/0-mercounts/$asm.ms$ms.histogram.info") or caExit("can't open '$wrk/0-mercounts/$asm.ms$ms.histogram.info' for reading: $!", undef);
            while (<F>) {
                $numTotal    = $1   if (m/Found\s(\d+)\s+mers./);
                $numDistinct = $1   if (m/Found\s(\d+)\s+distinct\smers./);
                $numUnique   = $1   if (m/Found\s(\d+)\s+unique\smers./);
                $largest     = $1   if (m/Largest\smercount\sis\s(\d+)/);
            }
            close(F);

            if (-e "$wrk/0-mercounts/$asm.ms$ms.histogram.sm.png") {
                push @$body, "<figure>\n";
                push @$body, "<a href='$wrk/0-mercounts/$asm.ms$ms.histogram.lg.png'><img src='$wrk/0-mercounts/$asm.ms$ms.histogram.sm.png'></a>\n";
                push @$body, "<figcaption>\n";
                push @$body, "Histogram for k=$ms with $numTotal mers, $numDistinct distinct mers and $numUnique single-copy mers.  Largest count is $largest.\n";
                push @$body, "</figcaption>\n";
                push @$body, "</figure>\n";
            }
        }

        elsif ((-e "$wrk/0-mercounts/$asm.ms$ms.ignore") && (-z "$wrk/0-mercounts/$asm.ms$ms.ignore")) {
            push @$body, "Threshold zdero.  No mers reported.\n";
        }

        elsif ((-e "$wrk/0-mercounts/$asm.ms$ms.fasta")  && (-z "$wrk/0-mercounts/$asm.ms$ms.fasta")) {
            push @$body, "Threshold zdero.  No mers reported.\n";
        }

        else {
            push @$body, "Using frequent mers in <some path>\n";
        }
    }
}



sub buildCorrectionHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Correction</h2>\n";
    push @$body, "\n";
    #buildGatekeeperHTML($wrk, $asm, $tag, $css, $body, $scripts);
    #  Analyzes the output fastq
}


sub buildCorrectedReadsHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Corrected Reads</h2>\n";
    push @$body, "\n";
    #buildGatekeeperHTML($wrk, $asm, $tag, $css, $body, $scripts);
    #  Analyzes the output fastq
}


sub buildTrimmingHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Trimming</h2>\n";
    push @$body, "\n";
    #buildGatekeeperHTML($wrk, $asm, $tag, $css, $body, $scripts);
    #  Analyzes the output fastq
}


sub buildTrimmedReadsHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Trimmed Reads</h2>\n";
    push @$body, "\n";
    #buildGatekeeperHTML($wrk, $asm, $tag, $css, $body, $scripts);
    #  Analyzes the output fastq
}


sub buildOverlapperHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Overlaps</h2>\n";
    push @$body, "\n";

    if (! -d "$wrk/$asm.ovlStore") {
        push @$body, "<p>Overlaps not computed.</p>\n";
        return;
    }

    if (! -e "$wrk/$asm.ovlStore.summary") {
        push @$body, "<p>No statistics available for store '$wrk/$asm.ovlStore'.</p>\n";
        return;
    }

    push @$body, "<table>\n";
    push @$body, "<tr><th>Category</th><th>Reads</th><th colspan='3'>Read Length</th><th colspan='3'>Feature Size or Coverage</th><th>Analysis</th></tr>\n";

    my ($category, $reads, $length, $lengthsd, $size, $sizesd, $analysis);

    open(F, "< $wrk/$asm.ovlStore.summary") or caExit("Failed to open overlap store statistics in '$wrk/$asm.ovlStore': $!", undef);
    $_ = <F>;
    $_ = <F>;
    while (<F>) {
        chomp;

        next if ($_ eq "");

        if      (m/(.*)\s+(\d+)\s+(\d+.\d+)\s+\+-\s+(\d+.\d+)\s+(\d+.\d+)\s+\+-\s+(\d+.\d+)\s+\((.*)\)$/) {
            $category = $1;
            $reads    = $2;
            $length   = $3;
            $lengthsd = $4;
            $size     = $5;
            $sizesd   = $6;
            $analysis = $7;
            push @$body, "<tr><td>$category</td><td>$reads</td><td align='right'>$length</td><td>&plusmn;</td><td align='left'>$lengthsd</td><td align='right'>$size</td><td>&plusmn;</td><td align='left'>$sizesd</td><td align='left'>$analysis</td></tr>\n";
            
        } elsif (m/(.*)\s+(\d+)\s+(\d+.\d+)\s+\+-\s+(\d+.\d+)\s+\((.*)\)$/) {
            $category = $1;
            $reads    = $2;
            $length   = $3;
            $lengthsd = $4;
            $size     = undef;
            $sizesd   = undef;
            $analysis = $5;
            push @$body, "<tr><td>$category</td><td>$reads</td><td align='right'>$length</td><td>&plusmn;</td><td align-'left'>$lengthsd</td><td></td><td></td><td></td><td align='left'>$analysis</td></tr>\n";

        } else {
            chomp;
            caExit("failed to parse line '$_' in file '$wrk/$asm.ovlStore.summary'", undef);
        }
    }
    close(F);

    push @$body, "</table>\n";
}


sub buildOverlapErrorCorrectionHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Overlap Error Adjustment</h2>\n";
    push @$body, "\n";
}


sub buildUnitiggerHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Unitigs</h2>\n";
    push @$body, "\n";
}


sub buildConsensusHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Consensus</h2>\n";
    push @$body, "\n";
}


sub buildOutputHTML ($$$$$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $css     = shift @_;  #  Array reference
    my $body    = shift @_;  #  Array reference
    my $scripts = shift @_;  #  Array reference

    push @$body, "<h2>Final Outputs</h2>\n";
    push @$body, "\n";
}


sub buildHTML ($$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $dir;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my @css;
    my @body;
    my @scripts;

    $dir = "correction"  if ($tag eq "cor");
    $dir = "trimming"    if ($tag eq "obt");
    $dir = "unitigging"  if ($tag eq "utg");

    $wrk = "$WRK/$dir";

    #  For correction runs
    if ($tag eq "cor") {
        push @body, "<h1>Correction</h1>\n";
        buildGatekeeperHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildMerylHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildOverlapperHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildCorrectionHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildCorrectedReadsHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
    }

    #  For trimming runs
    if ($tag eq "obt") {
        push @body, "<h1>Trimming</h1>\n";
        buildGatekeeperHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildMerylHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildOverlapperHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildTrimmingHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildTrimmedReadsHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
    }

    #  For assembly runs
    if ($tag eq "utg") {
        push @body, "<h1>Assembly</h1>\n";
        buildGatekeeperHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildMerylHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildOverlapperHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildOverlapErrorCorrectionHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildUnitiggerHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildConsensusHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
        buildOutputHTML($wrk, $asm, $tag, \@css, \@body, \@scripts);
    }


    #print STDERR "WRITING '$wrk/$asm-summary.html'\n";

    open(F, "> $WRK/$dir.html") or die "can't open '$WRK/$dir.html' for writing: $!\n";

    print F "<!DOCTYPE html>\n";
    print F "\n";
    print F "<html>\n";
    print F "\n";
    print F "<head>\n";
    print F "<title>canu analysis for assembly '$asm' in directory '$wrk'</title>\n";
    print F "<style type='text/css'>\n";
    print F "body       { font-family: Helvetica, Verdana, sans-serif; }\n";
    print F "h1, h2, h3 { color: #ee3e80; }\n";
    print F "p          { color: #665544; }\n";
    print F "th, td     { border: 1px solid #111111; padding: 2px 2px 2px 2px; }\n";
    print F "td:hover   { background-color: #e4e4e4; }\n";
    print F "th:hover   { background-color: #d4d4d4; }\n";
    print F "tr.details { visibility: collapse; }\n";
    print F @css;
    print F "</style>\n";
    print F "</head>\n";
    print F "\n";
    print F "<body>\n";
    print F "\n";
    print F @body;
    print F "\n";
    print F "<script type='text/javascript'>\n";
    print F "var toggleTable = function() {\n";
    print F "  var table = this.closest('table');\n";
    print F "  var elts  = table.querySelectorAll('.details');\n";
    print F "\n";
    print F "  for (var i=0; i<elts.length; i++) {\n";
    print F "    if (!elts[i].enabled) {\n";
    print F "      elts[i].enabled = true;\n";
    print F "      elts[i].style.visibility = 'visible';\n";
    print F "    } else {\n";
    print F "      elts[i].enabled = false;\n";
    print F "      elts[i].style.visibility = 'collapse';\n";
    print F "    }\n";
    print F "  }\n";
    print F "}\n";
    print F @scripts;
    print F "</script>\n";
    print F "\n";
    print F "</body>\n";
    print F "\n";
    print F "</html>\n";

    close(F);
}

