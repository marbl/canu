#!/usr/bin/perl

#  Runs the full assembly-to-assembly mapping comparison.
#
#  Takes two ATAC-format mapping files and:
#    shifts 1bp gaps to the 3' end
#    compute the assembly annotation
#
#  Generate statistics:
#    sum and histogram of the different annotations
#    sum and histogram of the disagreement 
#
#    number of matches in each
#    number of scaffolds mapped to multiple chromosomes
#
#    Nx (0 <= x <= 100)
#
#    histogram of match length
#    histogram of run length
#    histogram of clump length
#

my $bin             = "/bioinfo/assembly/walenz/hummap2/src/genomics/atac-driver";
my $overlap         = "$bin/alignOverlap/overlap";
my $gapShifter      = "$bin/gapShifter/gapShifter";
my $mismatchCounter = "$bin/mismatchCounter/mismatchCounter";

if (scalar(@ARGV) != 2) {
    print STDERR "usage: $0 <amap-mapping> <atac-mapping>\n";
    exit(1);
}

my $amap = shift @ARGV;
my $atac = shift @ARGV;
my $dir  = "COMPARE";
my $tmp  = "COMPARE";  #  /tmp on assembly-a is pathetic

#  If we're given relative paths, make them absolute
#
my $pwd = `pwd`;  chomp $pwd;

$amap = "$pwd/$amap" if ($amap !~ m!^/!);
$atac = "$pwd/$atac" if ($atac !~ m!^/!);


#  We run stuff in the $dir directory, set that up
#
system("mkdir $dir") if (! -d "$dir");
system("ln -s $amap $dir/amap.atac") if (! -e "$dir/amap.atac");
system("ln -s $atac $dir/atac.atac") if (! -e "$dir/atac.atac");


if (! -e "$dir/amap.shift.atac") {
    print STDERR "GAP SHIFTER NEEDS assembly-a, so that the whole genome can be loaded!\n";
    system("$gapShifter -g 1 < $dir/amap.atac > $dir/amap.shift.atac") and die "Failed to shift amap.\n";
}

if (! -e "$dir/atac.shift.atac") {
    print STDERR "GAP SHIFTER NEEDS assembly-a, so that the whole genome can be loaded!\n";
    system("$gapShifter -g 1 < $dir/atac.atac > $dir/atac.shift.atac") and die "Failed to shift atac.\n";
}


#  Run overlap
#
if (! -e "$dir/overlap.map1annotation") {
    system("$overlap $dir/amap.shift.atac $dir/atac.shift.atac $dir/overlap") and die "Failed to overlap.\n";
}


#  Reads the annotation output of overlap.C and adds percent identity
#  for each match.
#
#  This is done by first converting the overlap.C output into atac
#  format (one file for each mapping), running mismatchCounter on each
#  of those files, then merging the results together.
#
if (! -e "overlap.map1annotation.identity") {
    if ((-e "$tmp/c.otoi") && (-e "$tmp/d.otoi")) {
        print STDERR "Using $tmp/c.otoi and $tmp/d.otoi\n";
    } else {
        overlapToAtac("$dir/overlap.map1annotation", "$tmp/a.otoi", "$tmp/b.otoi");

        print STDERR "Counting mismatches.\n";

        system("$mismatchCounter < $tmp/a.otoi > $tmp/c.otoi") and die "Failed mismatchCounter on A.\n";
        system("$mismatchCounter < $tmp/b.otoi > $tmp/d.otoi") and die "Failed mismatchCounter on B.\n";
    }

    open(A, "< $tmp/c.otoi") or die "Failed to open $tmp/c.otoi\n";
    open(B, "< $tmp/d.otoi") or die "Failed to open $tmp/d.otoi\n";
    open(M, "< $dir/overlap.map1annotation");
    open(O, "> $dir/overlap.map1annotation.identity");

    my $a;
    my $b;
    my $m;

    print STDERR "Merging results.\n";

    while (!eof(A) && !eof(B) && !eof(M)) {
        $m = <M>;  chomp $m;

        #  Skip any ATAC headers in A and B
        do { $a = <A>;  chomp $a; } while ($a =~ m/^\//);
        do { $b = <B>;  chomp $b; } while ($b =~ m/^\//);

        my @av = split '\s+', $a;
        my @bv = split '\s+', $a;

        my ($aid, $abeg, $alen, $aori, $amis, $aident) = (undef, undef, undef, undef, undef, "0.0");
        my ($bid, $bbeg, $blen, $bori, $bmis, $bident) = (undef, undef, undef, undef, undef, "0.0");

        if ($a =~ m/HUREF:(\d+)\s+(\d+)\s+(\d+)\s+(-*\d+)\s+>\s+\/mismatches=(\d+)\s+\/identity=(\d+\.\d+)/) {
            $aid    = $1;
            $abeg   = $2;
            $alen   = $3;
            $aori   = $4;
            $amis   = $5;
            $aident = $6;
        } elsif ($a !~ m/^M\sm\s/) {
            print "Anope $a\n";
        }

        if ($b =~ m/HUREF:(\d+)\s+(\d+)\s+(\d+)\s+(-*\d+)\s+>\s+\/mismatches=(\d+)\s+\/identity=(\d+\.\d+)/) {
            $bid    = $1;
            $bbeg   = $2;
            $blen   = $3;
            $bori   = $4;
            $bmis   = $5;
            $bident = $6;
        } elsif ($b !~ m/^M\sm\s/) {
            print "Bnope $b\n";
        }

        $aident = substr("       $aident", -7);
        $bident = substr("       $bident", -7);

        if ($m =~ m/^(.*\]\s+\d+\s+\(.*\))\s+(\d+\s+\(.*\))/) {
            print O "$1 $aident $2 $bident\n";
        } else {
            print "Mnope $m\n";
            exit(1);
        }
    }
}




#  Reads the annotation output of overlap.C and writes two atac format
#  files, one for map1, one for map2.  There is a 1-1 map between
#  lines, unmapped B35 (either because B35 was unmapped by both
#  mappings, or unmapped by the other mapping) are noted in 'm'
#  matches ("M m UID . B35LC:xxxxx beg len 1")
#
sub overlapToAtac {
    my $infile  = shift @_;
    my $outfile1 = shift @_;
    my $outfile2 = shift @_;

    open(I,  "< $infile")  or die "Can't open '$infile' for reading.\n";
    open(O1, "> $outfile1") or die "Can't open '$outfile1' for writing.\n";
    open(O2, "> $outfile2") or die "Can't open '$outfile2' for writing.\n";

    print STDERR "Converting $infile -> $outfile1 and $outfile2\n";

    print O1 "/assemblyFile1=MERYL/B35LC.fasta\n";
    print O1 "/assemblyFile2=MERYL/HUREF2.fasta\n";
    print O1 "/assemblyId1=B35LC\n";
    print O1 "/assemblyId2=HUREF2\n";

    print O2 "/assemblyFile1=MERYL/B35LC.fasta\n";
    print O2 "/assemblyFile2=MERYL/HUREF2.fasta\n";
    print O2 "/assemblyId1=B35LC\n";
    print O2 "/assemblyId2=HUREF2\n";

    my $id = 0;

    while (<I>) {
        if (m/^.\s+(\d+):(\d+)-(\d+)\[\s*\d+\]\s(\d+)\s\(\s*(\d+):\s*(\d+)-\s*(\d+)\)\s(\d+)\s\(\s*(\d+):\s*(\d+)-\s*(\d+)\)\s*$/) {
            my $id1 = $1;
            my $b1  = $2;
            my $e1  = $3;
            my $l1  = $e1 - $b1;

            my $mid2a = $4;
            my $id2a  = $5;
            my $b2a   = $6;
            my $e2a   = $7;
            my $l2a   = $e2a - $b2a;
            my $oria  = 1;

            my $mid2b = $8;
            my $id2b  = $9;
            my $b2b   = $10;
            my $e2b   = $11;
            my $l2b   = $e2b - $b2b;
            my $orib  = 1;

            $b1 =~ s/^0+//;
            $e1 =~ s/^0+//;

            $b1 = 0 if ($b1 == 0);  #  fix for blowing away all the zeros

            if ($e2a < $b2a) {
                ($b2a, $e2a) = ($e2a, $b2a);
                $l2a         = $e2a - $b2a;
                $oria        = -1
                }
            if ($e2b < $b2b) {
                ($b2b, $e2b) = ($e2b, $b2b);
                $l2b         = $e2b - $b2b;
                $orib        = -1
                }

            $mid2a =~ s/^0+//;
            $mid2b =~ s/^0+//;

            if ($e2a > 0) {
                print O1 "M u $id . B35LC:$id1 $b1 $l1 1 HUREF:$id2a $b2a $l2a $oria\n";
            } else {
                print O1 "M m $id . B35LC:$id1 $b1 $l1 1\n";
            }

            if ($e2b > 0) {
                print O2 "M u $id . B35LC:$id1 $b1 $l1 1 HUREF:$id2b $b2b $l2b $orib\n";
            } else {
                print O2 "M m $id . B35LC:$id1 $b1 $l1 1\n";
            }

            $id++;
        } else {
            #print "Nope.\n";
            #exit(1);
        }
    }

    close(O1);
    close(O2);
    close(I);
}
