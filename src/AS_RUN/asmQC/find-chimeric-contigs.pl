#!/usr/bin/perl

use strict;

#  Given an assembly, map the contigs to a reference and report
#  errors -- chimerism, collapses, crud on the ends, etc.

my $snapper = "wgs/kmer/snapper/snapper2";
my $extent  = "wgs/kmer/sim4dbutils/convertToExtent";
my $leaff   = "wgs/kmer/leaff/leaff";
my $sim4th  = "wgs/kmer/sim4db/sim4th";

my $asm = shift @ARGV;
my $ref = shift @ARGV;

die "usage: $0 asm/asm.fasta ref.fasta\n" if (!defined($asm) || ! -e $asm);
die "usage: $0 asm/asm.fasta ref.fasta\n" if (!defined($ref) || ! -e $ref);

my %contigStatus;
my %contigLength;
my %contigBases;

open(F, "$leaff -F $asm -W |") or die;
while (!eof(F)) {
    my $h = <F>; chomp $h;
    my $s = <F>; chomp $s;

    if ($h =~ m/^>(\w+)/) {
        $contigStatus{$1} = "unmapped";
        $contigLength{$1} = length($s);
        $contigBases{$1}  = $s;
    } else {
        die "Failed to extract name from '$h'\n";
    }
}
close(F);

my $f = 0;
my $l = 0;
my $u = 0;


print STDERR "Mapping for full length matches.\n";
open(F, "$snapper -queries $asm -genomic $ref -noaligns | $extent -exons |");
$_ = <F>;  #  header line
while (<F>) {
    my @v = split '\s+', $_;

    if ($v[3] > $v[4]) {
        ($v[4], $v[3]) = ($v[3], $v[4]);
    }

    if ($v[4] - $v[3] == $v[1]) {
        print $_;
        $contigStatus{$v[0]} = "fulllength";
        $f++;
        $l += $v[1];
    }
}
close(F);

foreach my $c (keys %contigStatus) {
    if ($contigStatus{$c} eq "unmapped") {
        $u++;
    }
}

print STDERR "$f full length, total length $l, $u still unmapped.\n";

#  Build the next inputs

open(Z, "> find-chimer-$$.fasta");
foreach my $id (keys %contigStatus) {
    if ($contigStatus{$id} ne "fulllength") {
        print Z ">$id\n$contigBases{$id}\n";
    }
}
close(Z);

#  look for partial matches, filter out any small crap

#cDNAid  cDNAlen exonNum begin   end     genomicid       begin   end     identity        coverage
#ctg830  121541  0       121428  0       ntfl01_1        2828272 2949704 100.000  0.000
#ctg794  288336  0       288336  71791   ntfl01_1        758977  975529  100.000  0.000

my @allPartials;
my @goodPartials;

open(F, "$snapper -queries find-chimer-$$.fasta -genomic $ref -noaligns -minmatchcoverage 0 -discardexonlength 24 | $extent -exons |");
$_ = <F>;  #  header line
my @allPartials = <F>;
close(F);

unlink "find-chimer-$$.fasta";
unlink "find-chimer-$$.fastaidx";

print STDERR "found ", scalar(@allPartials), " partial matches.\n";

foreach my $tm (@allPartials) {
    my @t = split '\s+', $tm;
    ($t[3], $t[4]) = ($t[4], $t[3])  if ($t[3] > $t[4]);

    my @gp = @goodPartials;
    undef @goodPartials;

    my $doSave = 1;  #  We want to save this test match

    foreach my $gm (@gp) {
        my @g = split '\s+', $gm;
        ($g[3], $g[4]) = ($g[4], $g[3])  if ($g[3] > $g[4]);

        if ($g[0] ne $t[0]) {
            #  testMatch isn't from the same contig.  Resave the goodMatch
            push @goodPartials, $gm;
        } elsif (($g[3] >= $t[3]) && ($g[4] <= $t[4])) {
            #  former good is contained in new test, save new (unless flagged below)
        } elsif (($t[3] >= $g[3]) && ($t[4] <= $g[4])) {
            #  new test is contained in former good.  Resave the goodMatch.
            push @goodPartials, $gm;
            $doSave = 0;
        } else {
            #  Dovetail overlapping matches
            push @goodPartials, $gm;
        }
    }

    if ($doSave) {
        push @goodPartials, $tm;
    }
}

print STDERR "filtered down to ", scalar(@goodPartials), " partial matches.\n";

#open(F, "| sort -k7n");
#foreach my $gm (@goodPartials) {
#    print F $gm;
#}
#close(F);

@goodPartials = sort { $a cmp $b } @goodPartials;

chomp @goodPartials;

my %numMatches;

my $ctgid;
my @all;
my @agm;

my $badend   = 0;
my $collapse = 0;
my $chimer   = 0;
my $gaps     = 0;

push @goodPartials, "nothing 0 0 0 0 0 0 0 0 0\n";

foreach my $goodpart (@goodPartials) {
    my @v = split '\s+', $goodpart;

    ($v[3], $v[4]) = ($v[4], $v[3])  if ($v[3] > $v[4]);

    if ($ctgid ne $v[0]) {
        if (defined($ctgid)) {

            if (scalar(@all) == 1) {
                my ($ab, $ae, $rb, $re, $ln, $gm) = split '\0', $all[0];

                if ((0 < $ab) && ($ab < 1000)) {
                    #  Missing crud at the start
                    print "\n";
                    print "$gm\n";
                    print "*** crud at start\n";
                    $badend++;
                }

                if (($ae < $ln) && ($ae > $ln - 1000)) {
                    #  Missing crud at the end
                    print "\n";
                    print "$gm\n";
                    print "*** crud at end\n";
                    $badend++;
                }
            }

            if (scalar(@all) > 1) {
                @all = sort { $a <=> $b } @all;

                my @ab;
                my @ae;
                my @rb;
                my @re;
                my @gm;

                foreach my $a (@all) {
                    my ($ab, $ae, $rb, $re, $ln, $gm) = split '\0', $a;
                    push @ab, $ab;
                    push @ae, $ae;
                    push @rb, $rb;
                    push @re, $re;
                    push @gm, $gm;
                }

                my $i = 0;
                my $j = 1;

                while ($j < scalar(@ab)) {
                    my $asmDiff  = $ab[$j] - $ae[$i];
                    my $refDiffF = $rb[$j] - $re[$i];  #  assumes forward
                    my $refDiffR = $rb[$i] - $re[$j];  #  assumes reverse

                    print "\n";
                    print "$gm[$i]\n";
                    print "$gm[$j]\n";

                    if ($asmDiff == 0) {
                        #  No hole in assembly, what did the ref do?
                        if ($rb[$i] < $rb[$j]) {
                            print "NEGATIVE refDiffF $refDiffF\n" if ($refDiffF < 0);
                            if ($refDiffF < 1000) {
                                print "*** collapsed (forward) $refDiffF.\n";
                                $collapse++;
                                goto explained;
                            } else {
                                print "*** chimer (forward) $refDiffF.\n";
                                $chimer++;
                                goto explained;
                            }
                        } else {
                            print "NEGATIVE refDiffR $refDiffR\n" if ($refDiffR < 0);
                            if ($refDiffR < 1000) {
                                print "*** collapsed (reversed) $refDiffR.\n";
                                $collapse++;
                                goto explained;
                            } else {
                                print "*** chimer (forward) $refDiffR.\n";
                                $chimer++;
                                goto explained;
                            }
                        }
                    }

                    if ($asmDiff < 0) {
                        print "*** EXPANDED?\n";
                    }


                    #  else $asmDiff > 0
                    if ($rb[$i] < $rb[$j]) {
                        #  forward
                        print "*** gaps in both asm=$asmDiff ref=$refDiffF\n";
                        $gaps++;
                        goto explained;
                    } else {
                        #  reverse
                        print "*** gaps in both asm=$asmDiff ref=$refDiffR\n";
                        $gaps++;
                        goto explained;
                    }


                    print "*** unknown.\n";
                  explained:
                    $i++;
                    $j++;
                }
            }
        }

        $ctgid = $v[0];

        undef @all;
    }

    push @all, "$v[3]\0$v[4]\0$v[6]\0$v[7]\0$v[1]\0$goodpart";
}

print "badend=$badend collapse=$collapse chimer=$chimer gaps=$gaps\n";
