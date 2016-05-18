#!/usr/local/bin/perl

use strict;

my @dateStrings = ( "???", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" );

#  If set, rename original files to name.ORIG, rewrite files with updated copyright text.
#  If not, create new name.MODIFIED files with updated copyright text.
#
my $doForReal = 1;


#
#  The change data 'addCopyrights.dat' contains lines of two types:
#
#    A <file> <date><author>  -- a committed change to this file
#    D <file> <oldfile>       -- 'file' used to be called 'oldfile'
#
#  The 'D' lines DO NOT map existing A lines.  They just emit 'this file derived from' lines in the
#  copyright text.
#  
#  When a file is added to the repo, nothing special needs to be done; an 'A' line
#  is emitted for the new file.  If a file comes about because of a copy (of an existing file),
#  a 'D' line needs to be added.
#
#  If a file is renamed, a 'D' line needs to be added, and all previous mentions
#  of the original name need to be updated to the new name.
#


sub toList (@) {
    my @all = sort { $a <=> $b } @_;
    my $ret;

    $ret = substr($all[0], 0, 4);
    shift @all;

    foreach my $a (@all) {
        $a = substr($a, 0, 4);

        if ($ret =~ /^(\d+)$/) {
            if ($1 == $a) {
            } elsif ($1 + 1 == $a) {
                $ret = "$1-$a";
            } else {
                $ret = "$1,$a";
            }
        }

        if ($ret =~ /^(.*)-(\d+)$/) {
            if ($2 == $a) {
            } elsif ($2 + 1 == $a) {
                $ret = "$1-$a";
            } else {
                $ret = "$1-$2,$a";
            }
        }

        if ($ret =~ /^(.*),(\d+)$/) {
            if ($2 == $a) {
            } elsif ($2 + 1 == $a) {
                $ret = "$1,$2-$a";
            } else {
                $ret = "$1,$2,$a";
            }
        }
    }

    return($ret);
}




sub splitAC ($@) {
    my $cc = shift @_;
    my @AC = @_;
    my @AClist;

    my %dates;

    foreach my $ac (@AC) {
        if ($ac =~ m/^(....)(\d\d\d\d\d\d\d\d)(.*)$/) {
            $dates{"$1$3"} .= "$2\n";
        } else {
            die "$ac failed\n";
        }
    }

    foreach my $ac (keys %dates) {
        my @dates = split '\n', $dates{$ac};

        @dates = sort { $a <=> $b } @dates;

        my $years = toList(@dates);

        my $ord = $dates[0];
        my $bgn = $dates[0];
        my $end = $dates[ scalar(@dates)-1 ];

        if ($bgn =~ m/^(\d\d\d\d)(\d\d)(\d\d)$/) {
            $bgn = "$1-$dateStrings[$2]-$3";
        } else {
            die "bgn date $bgn\n";
        }

        if ($end =~ m/^(\d\d\d\d)(\d\d)(\d\d)$/) {
            $end = "$1-$dateStrings[$2]-$3";
        } else {
            die "bgn date $end\n";
        }

        my $org;
        my $nam;
        if ($ac =~ m/^(....)(.*)$/) {
            $org = $1;
            $nam = $2;
        } else {
            die "$ac match\n";
        }

        my $dates = "from $bgn to $end";
        if ($bgn eq $end) {
            $dates = "on $bgn";
        }

        if ($org eq "nihh") {
            $dates = "beginning on $bgn";
        }

        my $str;

        if      ($org eq "craa") {
            $str .= " $cc    $nam $dates\n";
            $str .= " $cc      are Copyright $years Applera Corporation, and\n";
            $str .= " $cc      are subject to the GNU General Public License version 2\n";
            $str .= " $cc\n";

        } elsif ($org eq "tigr") {
            $str .= " $cc    $nam $dates\n";
            $str .= " $cc      are Copyright $years The Institute for Genomics Research, and\n";
            $str .= " $cc      are subject to the GNU General Public License version 2\n";
            $str .= " $cc\n";

        } elsif ($org eq "jcvi") {
            $str .= " $cc    $nam $dates\n";
            $str .= " $cc      are Copyright $years J. Craig Venter Institute, and\n";
            $str .= " $cc      are subject to the GNU General Public License version 2\n";
            $str .= " $cc\n";

        } elsif ($org eq "bnbi") {
            $str .= " $cc    $nam $dates\n";
            $str .= " $cc      are Copyright $years Battelle National Biodefense Institute, and\n";
            $str .= " $cc      are subject to the BSD 3-Clause License\n";
            $str .= " $cc\n";

        } elsif ($org eq "nihh") {
            $str .= " $cc    $nam $dates\n";
            $str .= " $cc      are a 'United States Government Work', and\n";
            $str .= " $cc      are released in the public domain\n";
            $str .= " $cc\n";

        } elsif ($org eq "none") {
            $str .= " $cc    $nam $dates\n";
            $str .= " $cc      are Copyright $years $nam, and\n";
            $str .= " $cc      are subject to the GNU General Public License version 2\n";
            $str .= " $cc\n";

        } else {
            die "$ac org\n";
        }

        push @AClist, "$ord\0$str";
    }

    @AClist = sort { $a <=> $b } @AClist;

    foreach my $a (@AClist) {
        (undef, $a) = split '\0', $a;
    }

    return(@AClist);
}







#  Load the previously generated change data

my %authcopy;
my %derived;

{
    open(F, "< addCopyrights.dat");
    while (<F>) {
        chomp;

        if (m/^A\s+(\S+)\s+(.*)$/) {
            $authcopy{$1} .= "$2\n";

        } elsif (m/^D\s+(\S+)\s+(\S+)$/) {
            $derived{$1} .= "$2\n";

        } else {
            die "invalid addCopyrights.dat line '$_'\n";
        }
    }
    close(F);
}




#  Process each file.

open(FIN, "find kmer src -type f -print |") or die "Failed to launch 'find'\n";
#open(OUT, "> addCopyrights.dat.new") or die "Failed to open 'addCopyrights.dat.new' for writing: $!\n";

while (<FIN>) {
    chomp;

    my $file = $_;

    $file = $1  if ($_ =~ m/^\.\/(.*)$/);  #  Remove leading ./ added by find.

    my @lines;


    next if ($file =~ m/\.mk$/);
    next if ($file =~ m/Makefile/);

    next if ($file =~ m/\.sh$/);
    next if ($file =~ m/\.py$/);

    next if ($file =~ m/\.jar$/);
    next if ($file =~ m/\.tar$/);
    next if ($file =~ m/\.bin$/);  #  falcon_sense

    next if ($file =~ m/\.gz$/);

    next if ($file =~ m/\.json$/);
    next if ($file =~ m/\.json.README$/);
    next if ($file =~ m/\.css$/);

    next if ($file =~ m/\.fasta$/);  #  meryl test

    next if ($file =~ m/\.dat$/);  #  src/overlapInCore/liboverlap/prefixEditDistance-matchLimitData/prefixEditDistance-matchLimit-*.dat

    next if ($file =~ m/md5/);
    next if ($file =~ m/mt19937ar/);

    next if ($file =~ m/\.jpg$/);
    next if ($file =~ m/README/);
    next if ($file =~ m/\.dat$/);

    next if ($file =~ m/libboost/);

    my $cb = "/";
    my $cc = "*";
    my $ce = "/";

    if ($file =~ m/\.p[lm]$/) {
        $cb = "#";
        $cc = "#";
        $ce = "#";
    }

    my $iskmer   = 0;

    $iskmer = 1    if ($file =~ m/^kmer/);
    $iskmer = 1    if ($file =~ m/meryl/);

    $iskmer = 1    if ($file =~ m/bitEncodings/);
    $iskmer = 1    if ($file =~ m/bitOperations/);
    $iskmer = 1    if ($file =~ m/bitPackedArray/);
    $iskmer = 1    if ($file =~ m/bitPackedFile/);
    $iskmer = 1    if ($file =~ m/bitPacking/);
    $iskmer = 1    if ($file =~ m/decodeBooleanString/);
    $iskmer = 1    if ($file =~ m/dnaAlphabets/);
    $iskmer = 1    if ($file =~ m/intervalList/);
    $iskmer = 1    if ($file =~ m/kMer/);
    $iskmer = 1    if ($file =~ m/kMerHuge/);
    $iskmer = 1    if ($file =~ m/kMerTiny/);
    $iskmer = 1    if ($file =~ m/memoryMappedFile/);
    $iskmer = 1    if ($file =~ m/memoryMappedFileTest/);
    $iskmer = 1    if ($file =~ m/readBuffer/);
    $iskmer = 1    if ($file =~ m/speedCounter/);
    $iskmer = 1    if ($file =~ m/splitToWords/);
    $iskmer = 1    if ($file =~ m/sweatShop/);
    $iskmer = 1    if ($file =~ m/testHashTable/);
    $iskmer = 1    if ($file =~ m/testRand/);
    $iskmer = 1    if ($file =~ m/testVar/);
    $iskmer = 1    if ($file =~ m/timeAndSize/);

    if (($file !~ m/\.[CHch]$/) && ($file !~ m/\.p[lm]/)) {
        print STDERR "Won't process:      '$file'\n";
        next;
    }

    my @AC     = split '\n', $authcopy{$file};

    #foreach my $ac (@AC) {
    #    print OUT "A\t$file\t$ac\n";
    #}

    my @AClist = splitAC($cc, @AC);

    my %DElist;
    my @DElist = split '\n', $derived{$file};

    foreach my $d (@DElist) {
        next if ($d eq "");
        next if (lc $d eq lc $file);
        $DElist{$d}++;
    }

    undef @DElist;

    if (scalar(keys %DElist) > 0) {
        foreach my $d (keys %DElist) {
            push @DElist, " $cc    $d\n";
        }

        @DElist = sort @DElist;

        #foreach my $d (@DElist) {
        #    if ($d =~ m/^\s*\S\s*(\S+)$/) {
        #        print OUT "D\t$file\t$1\n";
        #    } else {
        #        die "Failed to match DElist line '$d'\n";
        #    }
        #}

        unshift @DElist, " $cc\n";
        unshift @DElist, " $cc  This file is derived from:\n";

        push    @DElist, " $cc\n";
    }

    if ($file =~ m/\.pl$/) {
        push @lines, "#!/usr/bin/env perl\n";
    }
    

    push @lines, "\n";
    push @lines, "$cb" . $cc x 78 . "\n"; 
    push @lines, " $cc\n";
    push @lines, " $cc  This file is part of canu, a software program that assembles whole-genome\n";
    push @lines, " $cc  sequencing reads into contigs.\n";
    push @lines, " $cc\n";
    push @lines, " $cc  This software is based on:\n";
    #push @lines, " $cc    RELEASE_1-3_2004-03-17 of the 'Celera Assembler' (http://wgs-assembler.sourceforge.net)\n";
    push @lines, " $cc    'Celera Assembler' (http://wgs-assembler.sourceforge.net)\n";
    push @lines, " $cc    the 'kmer package' (http://kmer.sourceforge.net)\n";
    push @lines, " $cc  both originally distributed by Applera Corporation under the GNU General\n";
    push @lines, " $cc  Public License, version 2.\n";
    push @lines, " $cc\n";
    push @lines, " $cc  Canu branched from Celera Assembler at its revision 4587.\n";
    push @lines, " $cc  Canu branched from the kmer project at its revision 1994.\n";
    push @lines, " $cc\n";
    push @lines, @DElist;
    push @lines, " $cc  Modifications by:\n";
    push @lines, " $cc\n";
    push @lines, @AClist;
    push @lines, " $cc  File 'README.licenses' in the root directory of this distribution contains\n";
    push @lines, " $cc  full conditions and disclaimers for each license.\n";
    push @lines, " $cc$ce\n";
    push @lines, "\n";

    my $start = 1;  #  To skip comment lines at the start of the file (the previous copyright block).

    open(F, "< $file") or die "Failed to open '$file' for reading: $!\n";
    while (<F>) {
        s/\s+$//;  #  Remove trailing spaces because they bug me.

        #  If a single "/*" at the start of the line, assume this is NOT an old copyright block,
        #  but a copyright block (or just a comment) from some third party code.

        if ($_ eq "\/\*") {
            print STDERR "Foreign code found: '$file'\n";
            $start = 0;
        }

        #  If not at the start, add the line.

        if ($start == 0) {
            push @lines, "$_\n";
            next;
        }

        #  Else, we're at the start; if blank or a comment, skip it.

        if (($_ eq "") ||            #  Blank lines
            ($_ =~ m/^[\/\s]\*/) ||  #  C-style comment (the old copyright block)
            ($_ =~ m/^\s*##/) ||     #  Perl comment, at least two #'s
            ($_ =~ m/^\s*#$/) ||     #  Perl comment, exactly one #
            ($_ =~ m/^\s*#\s/) ||    #  Perl comment, a single # followed by a space (so we don't catch #! lines)
            ($_ =~ m/^\s*#\!/)) {    #  #! line.  I guess we do want to skip them now
            next;
        }

        #  Else, add the line, and declare that we're no longer at the start.

        push @lines, "$_\n";

        $start = 0;
    }
    close(F);

    if ($doForReal) {
        rename "$file", "$file.ORIG";

        open(F, "> $file") or die "Failed to open '$file' for writing: $!\n";
        print F @lines;
        close(F);
    } else {
        open(F, "> $file.MODIFIED") or die "Failed to open '$file.MODIFIED' for writing: $!\n";
        print F @lines;
        close(F);
    }
}

close(FIN);
