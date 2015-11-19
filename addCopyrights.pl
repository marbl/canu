#!/usr/local/bin/perl

use strict;

my @dateStrings = ( "???", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" );

#
#  Lots of dead code (like those lists of names) from the first copyright add run.
#
#  This script should Just Run, if the last run was finished properly.  In particular, 'r6844' must be updated.
#
#  After a successful run, before committing:
#
#    update the log loop stopping condition to the revision you'll be committing.  Search for 'r6844'.
#    rename 'addCopyrights.dat.new' to 'addCopyrights.dat'.
#    CHECK YOUR WORK!  Diff a bunch of the changes.  Check new files.  Build it.  Does canu still run?
#    TEST COMMIT - commit this script, check that the commit shows up correctly.
#
#  Notes:
#
#    Most of the SVN log processing is useless now.  Github doesn't seem to report renames.
#

sub getName ($) {
    my $a = $_[0];
    my $A;

    if      ($a eq "adelcher")        {  $A = "Art Delcher";
    } elsif ($a eq "art.delcher")     {  $A = "Art Delcher";
    } elsif ($a eq "ahalpern")        {  $A = "Aaron Halpern";
    } elsif ($a eq "aaron.halpern")   {  $A = "Aaron Halpern";
    } elsif ($a eq "andreyto")        {
    } elsif ($a eq "andrey.tovchigrechko") {
    } elsif ($a eq "florea")          {  $A = "Liliana Florea";
    } elsif ($a eq "liliana.florea")  {  $A = "Liliana Florea";
    } elsif ($a eq "cmobarry")        {  $A = "Clark Mobarry";
    } elsif ($a eq "clark.mobarry")   {  $A = "Clark Mobarry";
    } elsif ($a eq "walenz")          {  $A = "Brian P. Walenz";
    } elsif ($a eq "bri")             {  $A = "Brian P. Walenz";
    } elsif ($a eq "brianwalenz")     {  $A = "Brian P. Walenz";
    } elsif ($a eq "brian p. walenz") {  $A = "Brian P. Walenz";
    } elsif ($a eq "brian.p..walenz") {  $A = "Brian P. Walenz";
    } elsif ($a eq "catmandew")       {  $A = "Ian Dew";
    } elsif ($a eq "ian.dew")         {  $A = "Ian Dew";
    } elsif ($a eq "eliv")            {  $A = "Eli Venter";
    } elsif ($a eq "eli.venter")      {  $A = "Eli Venter";
    } elsif ($a eq "gdenisov")        {  $A = "Gennady Denisov";
    } elsif ($a eq "gennady.denisov") {  $A = "Gennady Denisov";
    } elsif ($a eq "gesims")          {  $A = "Gregory Sims";
    } elsif ($a eq "greg.sims")       {  $A = "Gregory Sims";
    } elsif ($a eq "granger_sutton")  {  $A = "Granger Sutton";
    } elsif ($a eq "granger.sutton")  {  $A = "Granger Sutton";
    } elsif ($a eq "jason_miller")    {  $A = "Jason Miller";
    } elsif ($a eq "jasonmiller9704") {  $A = "Jason Miller";
    } elsif ($a eq "jason.miller")    {  $A = "Jason Miller";
    } elsif ($a eq "kli1000")         {  $A = "Kelvin Li";
    } elsif ($a eq "kelvin.li")       {  $A = "Kelvin Li";
    } elsif ($a eq "mcschatz")        {  $A = "Michael Schatz";
    } elsif ($a eq "michael.schatz")  {  $A = "Michael Schatz";
    } elsif ($a eq "mhayton")         {
    } elsif ($a eq "matt.hayton")     {
    } elsif ($a eq "mkotelbajcvi")    {
    } elsif ($a eq "m.kolteba")       {
    } elsif ($a eq "moweis")          {
    } elsif ($a eq "m.oweis")         {
    } elsif ($a eq "edwardnj")        {  #  kmer build system
    } elsif ($a eq "nathan.edwards")  {  #  kmer build system
    } elsif ($a eq "root")            {  #  Really?  possibly me on os-x
    } elsif ($a eq "halldobv")        {
    } elsif ($a eq "bjarni.halldorson") {
    } elsif ($a eq "fasulodp")        {  #  kmer build system
    } elsif ($a eq "dan.fasulo")      {  #  kmer build system
    } elsif ($a eq "rbolanos")        {  #  kmer build system
    } elsif ($a eq "randall.bolanos") {  #  kmer build system
    } elsif ($a eq "ripper")          {  #  kmer build system
    } elsif ($a eq "ross.lippert")    {  #  kmer build system
    } elsif ($a eq "skoren")          {  $A = "Sergey Koren";
    } elsif ($a eq "sergey.koren")    {  $A = "Sergey Koren";
    } elsif ($a eq "vrainish")        {
    } elsif ($a eq "v.rainish")       {
    } elsif ($a eq "walenzb")         {  $A = "Brian P. Walenz";
    } else {
        print STDERR "Unknown a '$a'";
    }

    return($A);
}



sub getCopyright ($$$$) {
    my $a = shift @_;
    my $y = shift @_;
    my $m = shift @_;
    my $d = shift @_;
    my $C;

    if      (($a eq "catmandew") && ($y eq "2004")) {         #  This is the initial commit.
    } elsif (($a eq "ian.dew")   && ($y eq "2004")) {         #  This is the initial commit.

    } elsif  ($a eq "adelcher")        {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "art.delcher")     {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "ahalpern")        {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "aaron.halpern")   {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "florea")          {  $C = "none$y$m$d";  #  Earliest 2010-07-08, Latest 2011-11-16
    } elsif  ($a eq "liliana.florea")  {  $C = "none$y$m$d";  #  Earliest 2010-07-08, Latest 2011-11-16
    } elsif  ($a eq "cmobarry")        {  $C = "craa$y$m$d";
    } elsif  ($a eq "clark.mobarry")   {  $C = "craa$y$m$d";
    } elsif  ($a eq "catmandew")       {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "ian.dew")         {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "eliv")            {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "eli.venter")      {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "gdenisov")        {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "gennady.denisov") {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "gesims")          {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "greg.sims")       {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "granger_sutton")  {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "granger.sutton")  {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "jason_miller")    {  $C = "tigr$y$m$d";
    } elsif  ($a eq "jasonmiller9704") {  $C = "jcvi$y$m$d";

    } elsif  ($a eq "jason.miller")    {
        if ($y < 2006) {
            $C = "tigr$y$m$d";
        } else {
            $C = "jcvi$y$m$d";
        }


    } elsif  ($a eq "kli1000")         {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "kelvin.li")       {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "mcschatz")        {  $C = "tigr$y$m$d";
    } elsif  ($a eq "michael.schatz")  {  $C = "tigr$y$m$d";

    } elsif ($a eq "sergey.koren") {
        if ($y <  "2011") {
            $C = "jcvi$y$m$d";
        } elsif (($y < "2015") || (($y eq "2015") && ($m < 10))) {
            $C = "bnbi$y$m$d";
        } else {
            $C = "nihh$y$m$d";
        }

    } elsif (($a eq "brianwalenz") || ($a eq "walenz") || ($a eq "bri") || ($a eq "walenzb") || ($a eq "brian.p..walenz")) {
        if      (($y < 2004) || (($y == 2004) && ($m < 10) && ($d < 9))) {
            $C = "craa$y$m$d";
        } elsif ($y < 2005) {
            $C = "none$y$m$d";  #  Not employed, copyright me
        } elsif (($y < 2014) || (($y == 2014) && ($m < 7))) {
            $C = "jcvi$y$m$d";
        } elsif (($y < 2015) || (($y < 2016) && ($m < 10))) {
            $C = "bnbi$y$m$d";
        } else {
            $C = "nihh$y$m$d";
        }

    } elsif ($a eq "mkotelbajcvi") {
    } elsif ($a eq "m.kolteba") {
    } elsif ($a eq "moweis") {
    } elsif ($a eq "m.oweis") {
    } elsif ($a eq "andreyto") {
    } elsif ($a eq "andrey.tovchigrechko") {
    } elsif ($a eq "vrainish") {
    } elsif ($a eq "v.rainish") {
    } elsif ($a eq "mhayton") {
    } elsif ($a eq "matt.hayton") {
    } elsif ($a eq "edwardnj") {
    } elsif ($a eq "nathan.edwards") {
    } elsif ($a eq "root") {
    } elsif ($a eq "halldobv") {
    } elsif ($a eq "bjarni.halldorson") {
    } elsif ($a eq "fasulodp") {
    } elsif ($a eq "dan.fasulo") {
    } elsif ($a eq "rbolanos") {
    } elsif ($a eq "randall.bolanos") {
    } elsif ($a eq "ripper") {
    } elsif ($a eq "ross.lippert") {

    } else {
        print STDERR "unknown name $a\n";
    }

    return($C);
}



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






#  Generate logs

if (! -e "logs") {
    my $rev = 0;

    open(F, "svn info |");
    while (<F>) {
        if (m/Revision:\s+(\d+)$/) {
            $rev = $1;
        }
    }
    close(F);

    open(L, "> logs");

    #  The first copyright addition is at r6844
    #  The second copyright addition is at r7001

    for (my $x=$rev; $x >= 6844; $x--) {
        print STDERR "$x\r";

        open(I, "svn log -v -r $x |");
        while (<I>) {
            print L $_;
        }
        close(I);
    }

    close(L);
}




my %authcopy;
my %derived;

{
    my $r = undef;
    my $a = undef;
    my $y = undef;
    my $m = undef;
    my $d = undef;

    #  Special case filemap; files that were split manually.
    #  Mostly this is from memory.
    #
    #  filemap{oldfile} -> newfile

    my %filemap;


    #  Read the logs, create data

    my $skip;

    open(F, "< logs") or die "Failed to open 'logs' for reading: $!\n";
    while (<F>) {
        chomp;

        if (m/^------------------------------------------------------------------------$/) {
            undef $r;
            undef $a;
            undef $y;
            undef $m;
            undef $d;
            undef $skip;

        } elsif (m/^\s*$/) {
            $skip = 1;

        } elsif ($skip) {

        } elsif (m/^r(\d+)\s\|\s(.*)\s\|\s(\d\d\d\d)-(\d\d)-(\d\d)/) {
            $r = $1;
            $a = $2;
            $y = $3;
            $m = $4;
            $d = $5;

        } elsif (m/^Changed paths:$/) {

        } elsif ($r == 1) {
            #  Skip the initial commit, it breaks the next clause.

        } elsif (m/^\s\s\s.\s\/trunk\/(\S+)/) {
            my $fn = $1;

            #  Add a filemap if this is a svn controlled move

            if      (m/^\s\s\s.\s\/trunk\/(\S+)\s.*from\s\/trunk\/(\S+)\:\d+/) {
                #print STDERR "SET '$2' -> '$1' (for '$fn')\n";
                if (defined($filemap{$2})) {
                    $filemap{$2} .= "\0$1";
                } else {
                    $filemap{$2} = $1;
                }
            } elsif (m/\s\sA\s\/trunk\/(\S+)\s/) {
                #print STDERR "NEW FILE $f\n";
            }

            #  And recursively figure out the files this change applies to in the future

            my %files;

            $files{$fn} = 1;  #  So far, just the current file.

            {
                my $added = 1;

                while ($added) {
                    $added = 0;

                    #  For every file currently in the list
                    foreach my $fk (keys %files) {

                        #  If it has a map, remove the old file, and add the mapped files.
                        if (exists($filemap{$fk})) {
                            #delete $files{$fk};

                            my @map = split '\0', $filemap{$fk};

                            foreach my $fm (@map) {
                                if (!exists($files{$fm})) {
                                    #print STDERR "ADD FILE for new $fm <- $fk <- $fn\n";
                                    $files{$fm} = 1;
                                    $added = 1;  #  If wasn't there, keep recursing
                                }
                            }
                        }
                    }
                }
            }

            my $A = getName($a);
            my $C = getCopyright($a, $y, $m, $d);

            if ((defined($A) && (defined($C)))) {
                foreach my $f (keys %files) {
                    $authcopy{$f}   .= "$C$A\n";
                    $derived{$f}    .= "$fn\n";    #  Whatever file we end up with, it was derived from the one we're reading the log for.
                }
            }
        } else {
            print STDERR "WARN: $_\n";
        }
    }
    close(F);

    #  Load the previous data

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
open(OUT, "> addCopyrights.dat.new") or die "Failed to open 'addCopyrights.dat.new' for writing: $!\n";

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


    my $cb = "/";
    my $cc = "*";
    my $ce = "/";

    if ($file =~ m/\.p[lm]$/) {
        $cb = "#";
        $cc = "#";
        $ce = "#";
    }



    my $iskmer   = 0;
    my $isextern = 0;

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

    $isextern = 1  if ($file =~ m/md5/);
    $isextern = 1  if ($file =~ m/md5/);

    $isextern = 1  if ($file =~ m/mt19937ar/);
    $isextern = 1  if ($file =~ m/mt19937ar/);


    die "Can't process '$file'\n"  if (($file !~ m/\.[CHch]$/) && ($file !~ m/\.p[lm]/));

    my @AC     = split '\n', $authcopy{$file};

    foreach my $ac (@AC) {
        print OUT "A\t$file\t$ac\n";
    }

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

        foreach my $d (@DElist) {
            if ($d =~ m/^\s*\S\s*(\S+)$/) {
                print OUT "D\t$file\t$1\n";
            } else {
                die "Failed to match DElist line '$d'\n";
            }
        }

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
            ($_ =~ m/^\s*#\!/)) {      #  #! line.  I guess we do want to skip them now
            next;
        }

        #  Else, add the line, and declare that we're no longer at the start.

        push @lines, "$_\n";

        $start = 0;
    }
    close(F);

    #rename "$file", "$file.ORIG";

    open(F, "> $file") or die "Failed to open '$file' for writing: $!\n";
    print F @lines;
    close(F);

    #open(F, "> $file.MODIFIED") or die "Failed to open '$file.MODIFIED' for writing: $!\n";
    #print F @lines;
    #close(F);
}

close(FIN);
