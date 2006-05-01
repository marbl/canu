#!/usr/bin/perl

#  Convers the output of statsGenerator into a nice Excel spreadsheet.

use strict;
use lib "/home/bwalenz/linux/lib/perl5/site_perl/5.8.0";
use Spreadsheet::WriteExcel;
use Spreadsheet::WriteExcel::Big;

if (scalar(@ARGV) != 1) {
    die "usage: $0 stats-prefix\n";
}

my $prefix = shift @ARGV;

################################################################################
#
#  First, suck in the big ugly stdout from statsGenerator.
#
if (! -e "$prefix.out") {
    die "I looked for the stdout from statsGenerator in '$prefix.out', but didn't find it.\n";
}

my $workbook = Spreadsheet::WriteExcel::Big->new("$prefix.xls");
my $summary  = $workbook->add_worksheet("Summary");

my $format = $workbook->add_format();
$format->set_size(10);
$format->set_color('black');
$format->set_num_format(1);

my $formatFP = $workbook->add_format();
$format->set_size(10);
$format->set_color('black');

my $format_heading = $workbook->add_format();
$format_heading->set_size(10);
$format_heading->set_bold();
$format_heading->set_color('black');
$format_heading->set_num_format(1);

my $format_label = $workbook->add_format();
$format_heading->set_size(10);
$format_heading->set_bold();
$format_heading->set_color('black');
$format_heading->set_num_format(1);

my $format_comment = $workbook->add_format();
$format_heading->set_size(10);
$format_heading->set_bold();
$format_heading->set_color('black');
$format_heading->set_num_format(1);

$summary->set_column(0, 2, 20);
$summary->set_column(3, 3, 30);

my %stats;  #  scratch space

open(F, "< $prefix.out");
while (!eof(F)) {
    $_ = <F>;

    if (m/^\s*$/) {
        #  Nop;
    } if (m/^SEQUENCE$/) {
        $summary->write(1, 0, "Input Sequences", $format_heading);
        $summary->write(2, 0, "totalLength", $format_label);
        $summary->write(2, 3, "all letters, including N", $format_comment);
        $summary->write(3, 0, "totalLength", $format_label);
        $summary->write(3, 3, "ACGT only", $format_comment);

        $_ = <F>;
        if (m/totalLength\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+#\s+all\s+letters/) {
            $summary->write(0, 1, $1, $format_heading);
            $summary->write(0, 2, $3, $format_heading);

            #  remember which column is for which assembly
            $stats{$1} = 1;
            $stats{$3} = 2;

            #  and which assembly is in which column
            $stats{1} = $1;
            $stats{2} = $3;

            $summary->write(2, 1, "$2", $format);
            $summary->write(2, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }

        $_ = <F>;
        if (m/totalLength\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+#\s+ACGT/) {
            $summary->write(3, 1, "$2", $format);
            $summary->write(3, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }
    } if (m/^TANDEM REPEATS in (.*)$/) {
        my $asm = $1;

        $summary->write(5, 0, "Tandem Repeats", $format_heading);
        $summary->write(6, 0, "number", $format_label);
        $summary->write(7, 0, "totalLength", $format_label);
        $summary->write(8, 0, "coveredLength", $format_label);

        $_ = <F>;
        if (m/numberOfItems\s+(\d+)/) {
            $summary->write(6, $stats{$asm}, "$1", $format);
        } else {
            die "Parse error $_\n";
        }

        $_ = <F>;
        if (m/totalLength\s+(\d+)\s+#\s+sum\s+of\s+lengths/) {
            $summary->write(7, $stats{$asm}, "$1", $format);
        } else {
            die "Parse error $_\n";
        }

        $_ = <F>;
        if (m/coveredLength\s+(\d+)\s+#\s+sequence\s+covered/) {
            $summary->write(8, $stats{$asm}, "$1", $format);
        } else {
            die "Parse error $_\n";
        }
    } if (m/^MATCHES IN RUNS$/) {
        $summary->write(10, 0, "Matches in Runs", $format_heading);
        $summary->write(11, 0, "runMissingFull", $format_label);
        $summary->write(11, 3, "covered by a run, not by a match, including N", $format_comment);
        $summary->write(12, 0, "runMissingACGT", $format_label);
        $summary->write(12, 3, "covered by a run, not by a match, ACGT only", $format_comment);

        $_ = <F>;
        if (m/runMissingFull\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+#\s+sequence\s+in\s+run,\s+not\s+covered,\s+including\s+N/) {
            $summary->write(11, 1, "$2", $format);
            $summary->write(11, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }

        $_ = <F>;
        if (m/runMissingFull\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+#\s+sequence\s+in\s+run,\s+not\s+covered,\s+ACGT\s+only/) {
            $summary->write(12, 1, "$2", $format);
            $summary->write(12, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }
    } if ((m/^MATCHES$/) || (m/^RUNS$/)) {
        my $begin;
        my $chrcov;

        if (m/^MATCHES$/) {
            $begin = 14;
            $summary->write($begin, 0, "Matches", $format_heading);
            $chrcov  = $workbook->add_worksheet("Chr Cov Match");
        } else {
            $begin = 26;
            $summary->write($begin, 0, "Runs", $format_heading);
            $chrcov  = $workbook->add_worksheet("Chr Cov Run");
        }

        $chrcov->set_column(0, 0, 10);
        $chrcov->set_column(1, 6, 20);

        $summary->write($begin+1, 0, "number", $format_label);
        $summary->write($begin+2, 0, "totalLength", $format_label);
        $summary->write($begin+3, 1, "histogram", $format_label);
        $summary->write($begin+4, 2, "histogram", $format_label);

        $summary->write($begin+5, 0, "coveredLengthFull", $format_label);
        $summary->write($begin+6, 0, "coveredLengthACGT", $format_label);
        $summary->write($begin+7, 0, "coveredLengthNonACGT", $format_label);
        $summary->write($begin+8, 1, "histogram", $format_label);
        $summary->write($begin+9, 2, "histogram", $format_label);

        $_ = <F>;
        if (m/numberOfItems\s+(\d+)/) {
            $summary->write($begin+1, 1, "$1", $format);
        } else {
            die "Parse error $_\n";
        }
        
        $_ = <F>;
        if (m/matchLength\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s#\sSum\s+of\s+lengths/) {
            $summary->write($begin+2, 1, "$2", $format);
            $summary->write($begin+2, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }
        
        histogram($begin+3);
        histogram($begin+4);

        $_ = <F>;
        if (m/coveredLength\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s/) {
            $summary->write($begin+5, 1, "$2", $format);
            $summary->write($begin+5, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }
        
        $_ = <F>;
        if (m/coveredLength\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s/) {
            $summary->write($begin+6, 1, "$2", $format);
            $summary->write($begin+6, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }

        $_ = <F>;
        if (m/coveredLength\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s/) {
            $summary->write($begin+7, 1, "$2", $format);
            $summary->write($begin+7, 2, "$4", $format);
        } else {
            die "Parse error $_\n";
        }

        histogram($begin+8);
        histogram($begin+9);

        #  chromosome covered

        $chrcov->write(0, 1, "all sequence covered", $format_heading);
        $chrcov->write(0, 2, "all sequence length", $format_heading);
        $chrcov->write(0, 3, "percent covered", $format_heading);

        $chrcov->write(0, 4, "ACGT sequence covered", $format_heading);
        $chrcov->write(0, 5, "ACGT sequence length", $format_heading);
        $chrcov->write(0, 6, "percent covered", $format_heading);

        for (my $i=1; $i<23; $i++) {
            $chrcov->write($i, 0, "Chr$i", $format_heading);
        }
        $chrcov->write(23, 0, "ChrX", $format_heading);
        $chrcov->write(24, 0, "ChrY", $format_heading);
        $chrcov->write(25, 0, "ChrMT", $format_heading);
        $chrcov->write(26, 0, "ChrUn", $format_heading);

        $_ = <F>;
        while (m/chrCoveredLength\[\s*(\d+)\]\s+\S+\s+(\d+)\s+(\d+)\s+(\d+.\d+)%\s+(\d+)\s+(\d+)\s+(\d+.\d+)%\s+/) {
            my $i = $1 + 1;
            my $j = $1 + 2;

            $chrcov->write($i, 1, "$2", $format);
            $chrcov->write($i, 2, "$3", $format);
            $chrcov->write($i, 3, "=B$j / C$j * 100", $formatFP);
            $chrcov->write($i, 4, "$5", $format);
            $chrcov->write($i, 5, "$6", $format);
            $chrcov->write($i, 6, "=E$j / F$j * 100", $formatFP);

            $_ = <F>;
        }
    }
}
close(F);



sub histogram ($) {
    my $x = shift @_;

    $_ = <F>;
    if (m/histogram\s+\S+\s+(\d+)\s+items\s+(\d+.\d+)\s+average\s+(\d+.\d+)\s+std.dev./) {
        $summary->write($x, 3, "items, average, std.dev.", $format_comment);
        $summary->write($x, 4, "$1", $format);
        $summary->write($x, 5, "$2", $formatFP);
        $summary->write($x, 6, "$3", $formatFP);
    } else {
        die "Parse error $_\n";
    }
}


################################################################################
#
#  Nx
#
my $nx = $workbook->add_worksheet("Nx");

$nx->write(0, 1, "matches", $format_heading);
$nx->write(0, 2, "runs", $format_heading);

open(A, "< $prefix-matches.Nx");
open(B, "< $prefix-runs.Nx");
while (!eof(A)) {
    my $a = <A>;
    my $b = <B>;
    my ($ai, $an) = split '\s+', $a;
    my ($bi, $bn) = split '\s+', $b;
    die "Nx error: ai=$ai != bi=$bi\n" if ($ai != $bi);
    $nx->write($ai, 0, $ai);
    $nx->write($ai, 1, $an);
    $nx->write($ai, 2, $bn);
}
close(B);
close(A);


################################################################################
#
#  Histograms of lengths
#
#
my $sheet = $workbook->add_worksheet("Matches Histogram");
$sheet->set_column(0, 4, 25);
$sheet->write(0, 1, "$stats{1} length", $format_heading);
$sheet->write(0, 2, "$stats{2} length", $format_heading);
$sheet->write(0, 3, "$stats{1} covered N", $format_heading);
$sheet->write(0, 4, "$stats{2} covered N", $format_heading);

dumpHistogram($sheet,
              "$prefix-matches.AmatchLength.histogramdat",
              "$prefix-matches.BmatchLength.histogramdat",
              "$prefix-matches.AcoveredN.histogramdat",
              "$prefix-matches.BcoveredN.histogramdat");


my $sheet = $workbook->add_worksheet("Runs Histogram");
$sheet->set_column(0, 4, 25);
$sheet->write(0, 1, "$stats{1} length", $format_heading);
$sheet->write(0, 2, "$stats{2} length", $format_heading);
$sheet->write(0, 3, "$stats{1} covered N", $format_heading);
$sheet->write(0, 4, "$stats{2} covered N", $format_heading);

dumpHistogram($sheet,
              "stats-runs.AmatchLength.histogramdat",
              "stats-runs.BmatchLength.histogramdat",
              "stats-runs.AcoveredN.histogramdat",
              "stats-runs.BcoveredN.histogramdat");

my $sheet = $workbook->add_worksheet("Run Missing Histogram");
$sheet->set_column(0, 4, 25);
$sheet->write(0, 1, "$stats{1} full missing", $format_heading);
$sheet->write(0, 3, "$stats{2} full missing", $format_heading);
$sheet->write(0, 2, "$stats{1} ACGT missing", $format_heading);
$sheet->write(0, 4, "$stats{2} ACGT missing", $format_heading);

dumpHistogram($sheet,
              "stats.ARunMissingFull.histogramdat",
              "stats.BRunMissingFull.histogramdat",
              "stats.ARunMissingACGT.histogramdat",
              "stats.BRunMissingACGT.histogramdat");



sub dumpHistogram {
    my $sheet = shift @_;
    my @files = @_;
    my $col   = 1;
    my $idx   = 0;
    my @range;

    #  I can't seem to find any way of deleting a cell once it's
    #  written (opposed to simply clearing the cell).  We want
    #  to know what the maximum value in any histogram is, so
    #  we can stop writing after that point.
    #
    my $maxIdx = 0;

    foreach my $f (@files) {
        $idx = 0;

        open(F, "< $f") or die "Failed to open dumpHistogram1 '$f'\n";
        my @lines = <F>;
        close(F);

        #  Don't use the last line in the file - this is the number
        #  of things bigger than the max, we always report this.
        pop @lines;

        foreach my $l (@lines) {
            my ($r, $v) = split '\s+', $l;
            $maxIdx = $idx if ($v > 0) && ($maxIdx < $idx);
            $idx++;
        }
    }

    #  Read the range from the first file -- we'll check that all the other files
    #  use the same range.
    #
    $idx = 0;
    open(F, "< $files[0]") or die "Failed to open dumpHistogram2 '$files[0]'\n";
    while (<F>) {
        my ($r, $v) = split '\s+', $_;
        $range[$idx] = $r;
        $sheet->write($idx+1, 0, "$r", $format);
        $idx++;
        last if ($idx > $maxIdx);
    }
    my $lastVal;
    while (<F>) {
        my ($r, $v) = split '\s+', $_;
        $lastVal = $r;
    }
    $sheet->write($idx+1, 0, "$lastVal", $format);
    close(F);



    foreach my $f (@files) {
        $idx = 0;

        open(F, "< $f") or die "Failed to open dumpHistogram3 '$f'\n";
        while (<F>) {
            my ($r, $v) = split '\s+', $_;
            die "range error in file '$f' at idx $idx; $range[$idx] != $r\n" if ($range[$idx] != $r);
            $sheet->write($idx+1, $col, "$v", $format);
            $idx++;
            last if ($idx > $maxIdx);
        }
        my $lastVal;
        while (<F>) {
            my ($r, $v) = split '\s+', $_;
            $lastVal = $v;
        }
        $sheet->write($idx+1, $col, "$lastVal", $format);
        close(F);

        $col++;
    }
}



################################################################################
#
#  By chromosome histograms of lengths
#

my $sheet = $workbook->add_worksheet("Matches Chr Histogram");
$sheet->set_column(0, 32, 12);

for (my $i=1; $i<23; $i++) {
    $sheet->write(0, $i, "Chr$i ACGT", $format_heading);
}
$sheet->write(0, 23, "ChrX ACGT", $format_heading);
$sheet->write(0, 24, "ChrY ACGT", $format_heading);
$sheet->write(0, 25, "ChrMT ACGT", $format_heading);
$sheet->write(0, 26, "ChrUn ACGT", $format_heading);

dumpHistogram($sheet,
              "stats-matches.chr00acgt.histogramdat",
              "stats-matches.chr01acgt.histogramdat",
              "stats-matches.chr02acgt.histogramdat",
              "stats-matches.chr03acgt.histogramdat",
              "stats-matches.chr04acgt.histogramdat",
              "stats-matches.chr05acgt.histogramdat",
              "stats-matches.chr06acgt.histogramdat",
              "stats-matches.chr07acgt.histogramdat",
              "stats-matches.chr08acgt.histogramdat",
              "stats-matches.chr09acgt.histogramdat",
              "stats-matches.chr10acgt.histogramdat",
              "stats-matches.chr11acgt.histogramdat",
              "stats-matches.chr12acgt.histogramdat",
              "stats-matches.chr13acgt.histogramdat",
              "stats-matches.chr14acgt.histogramdat",
              "stats-matches.chr15acgt.histogramdat",
              "stats-matches.chr16acgt.histogramdat",
              "stats-matches.chr17acgt.histogramdat",
              "stats-matches.chr18acgt.histogramdat",
              "stats-matches.chr19acgt.histogramdat",
              "stats-matches.chr20acgt.histogramdat",
              "stats-matches.chr21acgt.histogramdat",
              "stats-matches.chr22acgt.histogramdat",
              "stats-matches.chr23acgt.histogramdat",
              "stats-matches.chr24acgt.histogramdat",
              "stats-matches.chr25acgt.histogramdat");


my $sheet = $workbook->add_worksheet("Runs Chr ACGT Histogram");
$sheet->set_column(0, 32, 12);

for (my $i=1; $i<23; $i++) {
    $sheet->write(0, $i, "Chr$i ACGT", $format_heading);
}
$sheet->write(0, 23, "ChrX ACGT", $format_heading);
$sheet->write(0, 24, "ChrY ACGT", $format_heading);
$sheet->write(0, 25, "ChrMT ACGT", $format_heading);
$sheet->write(0, 26, "ChrUn ACGT", $format_heading);

dumpHistogram($sheet,
              "stats-runs.chr00acgt.histogramdat",
              "stats-runs.chr01acgt.histogramdat",
              "stats-runs.chr02acgt.histogramdat",
              "stats-runs.chr03acgt.histogramdat",
              "stats-runs.chr04acgt.histogramdat",
              "stats-runs.chr05acgt.histogramdat",
              "stats-runs.chr06acgt.histogramdat",
              "stats-runs.chr07acgt.histogramdat",
              "stats-runs.chr08acgt.histogramdat",
              "stats-runs.chr09acgt.histogramdat",
              "stats-runs.chr10acgt.histogramdat",
              "stats-runs.chr11acgt.histogramdat",
              "stats-runs.chr12acgt.histogramdat",
              "stats-runs.chr13acgt.histogramdat",
              "stats-runs.chr14acgt.histogramdat",
              "stats-runs.chr15acgt.histogramdat",
              "stats-runs.chr16acgt.histogramdat",
              "stats-runs.chr17acgt.histogramdat",
              "stats-runs.chr18acgt.histogramdat",
              "stats-runs.chr19acgt.histogramdat",
              "stats-runs.chr20acgt.histogramdat",
              "stats-runs.chr21acgt.histogramdat",
              "stats-runs.chr22acgt.histogramdat",
              "stats-runs.chr23acgt.histogramdat",
              "stats-runs.chr24acgt.histogramdat",
              "stats-runs.chr25acgt.histogramdat");




my $sheet = $workbook->add_worksheet("Runs Chr Full Histogram");
$sheet->set_column(0, 32, 12);

for (my $i=1; $i<23; $i++) {
    $sheet->write(0, $i, "Chr$i Full", $format_heading);
}
$sheet->write(0, 23, "ChrX FULL", $format_heading);
$sheet->write(0, 24, "ChrY FULL", $format_heading);
$sheet->write(0, 25, "ChrMT FULL", $format_heading);
$sheet->write(0, 26, "ChrUn FULL", $format_heading);

dumpHistogram($sheet,
              "stats-runs.chr00full.histogramdat",
              "stats-runs.chr01full.histogramdat",
              "stats-runs.chr02full.histogramdat",
              "stats-runs.chr03full.histogramdat",
              "stats-runs.chr04full.histogramdat",
              "stats-runs.chr05full.histogramdat",
              "stats-runs.chr06full.histogramdat",
              "stats-runs.chr07full.histogramdat",
              "stats-runs.chr08full.histogramdat",
              "stats-runs.chr09full.histogramdat",
              "stats-runs.chr10full.histogramdat",
              "stats-runs.chr11full.histogramdat",
              "stats-runs.chr12full.histogramdat",
              "stats-runs.chr13full.histogramdat",
              "stats-runs.chr14full.histogramdat",
              "stats-runs.chr15full.histogramdat",
              "stats-runs.chr16full.histogramdat",
              "stats-runs.chr17full.histogramdat",
              "stats-runs.chr18full.histogramdat",
              "stats-runs.chr19full.histogramdat",
              "stats-runs.chr20full.histogramdat",
              "stats-runs.chr21full.histogramdat",
              "stats-runs.chr22full.histogramdat",
              "stats-runs.chr23full.histogramdat",
              "stats-runs.chr24full.histogramdat",
              "stats-runs.chr25full.histogramdat");


