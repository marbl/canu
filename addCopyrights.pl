#!/usr/local/bin/perl

use strict;


sub loadFile ($@) {
    my $file  = shift @_;
    my @lines = @_;

    my $start = 1;

    open(F, "< $file") or die "Failed to open '$file' for reading: $!\n";
    while (<F>) {
        s/\s+$//;

        if ($start == 0) {
            push @lines, "$_\n";
            next;
        }

        if (($_ eq "") || ($_ =~ m/^[\/\s]\*/)) {
            #print STDERR "$file -- $_\n";
            next;
        }

        push @lines, "$_\n";

        $start = 0;
    }
    close(F);

    return(@lines);
}


sub getName ($) {
    my $a = $_[0];
    my $A;

    if      ($a eq "adelcher") {
        $A = "Art Delcher";

    } elsif ($a eq "ahalpern") {
        $A = "Aaron Halpern";

    } elsif ($a eq "andreyto") {

    } elsif ($a eq "florea") {
        $A = "Liliana Florea";
    } elsif ($a eq "cmobarry") {
        $A = "Clark Mobarry";

    } elsif ($a eq "walenz") {
        $A = "Brian P. Walenz";

    } elsif ($a eq "bri") {
        $A = "Brian P. Walenz";

    } elsif ($a eq "brianwalenz") {
        $A = "Brian P. Walenz";

    } elsif ($a eq "catmandew") {
        $A = "Ian Dew";

    } elsif ($a eq "eliv") {
        $A = "Eli Venter";

    } elsif ($a eq "gdenisov") {
        $A = "Gennady Denisov";

    } elsif ($a eq "gesims") {
        $A = "Gregory Sims";

    } elsif ($a eq "granger_sutton") {
        $A = "Granger Sutton";

    } elsif ($a eq "jason_miller") {
        $A = "Jason Miller";

    } elsif ($a eq "jasonmiller9704") {
        $A = "Jason Miller";

    } elsif ($a eq "kli1000") {
        $A = "Kelvin Li";

    } elsif ($a eq "mcschatz") {
        $A = "Michael Schatz";

    } elsif ($a eq "mhayton") {

    } elsif ($a eq "mkotelbajcvi") {

    } elsif ($a eq "moweis") {

    } elsif ($a eq "edwardnj") {
        #  kmer build system

    } elsif ($a eq "root") {
        #  Really?  possibly me on os-x

    } elsif ($a eq "halldobv") {

    } elsif ($a eq "fasulodp") {
        #  kmer build system

    } elsif ($a eq "rbolanos") {
        #  kmer build system

    } elsif ($a eq "ripper") {
        #  kmer build system

    } elsif ($a eq "skoren") {
        $A = "Sergey Koren";

    } elsif ($a eq "vrainish") {

    } elsif ($a eq "walenzb") {
        $A = "Brian P. Walenz";

    } else {
        die "Unknown a '$a'";
    }

    return($A);
}



sub getCopyright ($$$$) {
    my $a = shift @_;
    my $y = shift @_;
    my $m = shift @_;
    my $d = shift @_;

    my $C;

    if      ($a eq "adelcher") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "ahalpern") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "florea") {
        #  Earliest 2010-07-08
        #  Latest   2011-11-16
        $C = "none$y$m$d";

    } elsif ($a eq "cmobarry") {
        #  Earliest 2010-07-08
        #  Latest   2011-11-16
        $C = "craa$y$m$d";

    } elsif (($a eq "brianwalenz") || ($a eq "walenz") || ($a eq "bri") || ($a eq "walenzb")) {
        if      (($y < 2004) || (($y == 2004) && ($m < 10) && ($d < 9))) {
            $C = "craa$y$m$d";
        } elsif ($y < 2005) {
            $C = "none$y$m$d";  #  Not employed, copyright me
        } elsif (($y < 2014) || (($y == 2014) && ($m < 7))) {
            $C = "jcvi$y$m$d";
        } else {
            $C = "bnbi$y$m$d";
        }

    } elsif (($a eq "catmandew") && ($y eq "2004")) {
        #  This is the initial commit.
        #$C = "craa$y$m$d";

    } elsif ($a eq "catmandew") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "eliv") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "gdenisov") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "gesims") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "granger_sutton") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "jason_miller") {
        $C = "tigr$y$m$d";

    } elsif ($a eq "jasonmiller9704") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "kli1000") {
        $C = "jcvi$y$m$d";

    } elsif ($a eq "mcschatz") {
        $C = "tigr$y$m$d";

    } elsif ($a eq "skoren") {
        if ($y < "2011") {
            $C = "jcvi$y$m$d";
        } else {
            $C = "bnbi$y$m$d";
        }

        #  Ignores
    } elsif ($a eq "mkotelbajcvi") {
    } elsif ($a eq "moweis") {
    } elsif ($a eq "andreyto") {
    } elsif ($a eq "vrainish") {
    } elsif ($a eq "mhayton") {
    } elsif ($a eq "edwardnj") {
    } elsif ($a eq "root") {
    } elsif ($a eq "halldobv") {
    } elsif ($a eq "fasulodp") {
    } elsif ($a eq "rbolanos") {
    } elsif ($a eq "ripper") {

    } else {
        die "unknown name $a\n";
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





sub splitA (@) {
    my @A = @_;
    my %A;
    my @Alist;

    foreach my $a (@A) {
        $A{$a}++;
    }
    foreach my $a (keys %A) {
        my $com = ($A{$a} == 1) ? "commit" : "commits";

        my $c = substr(" *   $a ($A{$a} $com)                                                                                ", 0, 98) . "*\n";
        push @Alist, $c;
    }

    return(@Alist);
}



sub splitC (@) {
    my @C = @_;

    my %CRAAcounts;
    my %TIGRcounts;
    my %JCVIcounts;
    my %BNBIcounts;

    foreach my $c (@C) {
        if ($c =~ m/^(....)(\d\d\d\d)\d\d\d\d$/) {
            if      ($1 eq "craa") {
                $CRAAcounts{$2}++;

            } elsif ($1 eq "tigr") {
                $TIGRcounts{$2}++;

            } elsif ($1 eq "jcvi") {
                $JCVIcounts{$2}++;

            } elsif ($1 eq "bnbi") {
                $BNBIcounts{$2}++;

            } else {
                die "$c name\n";
            }
        } else {
            die "$c match\n";
        }

        #$C{$c}++;
    }

    my @Clist;

    if (scalar(keys %CRAAcounts) > 0) {
        my $d = toList(keys %CRAAcounts);
        my $c = " * Copyright (C) $d, Applera Corporation.";
        $c = substr($c . "                                                                                ", 0, 98) . "*\n";
        push @Clist, $c;
    }

    if (scalar(keys %TIGRcounts) > 0) {
        my $d = toList(keys %TIGRcounts);
        my $c = " * Copyright (C) $d, The Institute for Genomics Research.";
        $c = substr($c . "                                                                                ", 0, 98) . "*\n";
        push @Clist, $c;
    }

    if (scalar(keys %JCVIcounts) > 0) {
        my $d = toList(keys %JCVIcounts);
        my $c = " * Copyright (C) $d, J. Craig Venter Institute.";
        $c = substr($c . "                                                                                ", 0, 98) . "*\n";
        push @Clist, $c;
    }

    if (scalar(keys %BNBIcounts) > 0) {
        my $d = toList(keys %BNBIcounts);
        my $c = " * Copyright (C) $d, Battelle National Biodefense Institute.";
        $c = substr($c . "                                                                                ", 0, 98) . "*\n";
        push @Clist, $c;
    }

    return(@Clist);
}


#  Returns
#    Modifications by BPW from date-date are Copyright 2014-2019 J. Craig Venter Institute, and are covered under the General Public License
#
#  needs to build list of {name}{place} with the dates mods were made, then for each {name}{place} to find the min/max.
#

my @dateStrings = ( "???", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" );

sub splitAC (@) {
    my @AC = @_;

    my %dates;

    foreach my $ac (@AC) {
        if ($ac =~ m/^(....)(\d\d\d\d\d\d\d\d)(.*)$/) {
            $dates{"$1$3"} .= "$2\n";
        } else {
            die "$ac failed\n";
        }
    }

    #  Figure out the longest year string, so we can keep columns

    #my $yearlen = 0;
    #
    #foreach my $ac (keys %dates) {
    #    my @dates = split '\n', $dates{$ac};
    #
    #    @dates = sort { $a <=> $b } @dates;
    #
    #    my $years = toList(@dates);
    #
    #    if ($yearlen < length($years)) {
    #        $yearlen = length($years);
    #    }
    #}


    my @AClist;

    foreach my $ac (keys %dates) {
        my @dates = split '\n', $dates{$ac};

        @dates = sort { $a <=> $b } @dates;

        #my $years = substr(toList(@dates) . " " x $yearlen, 0, $yearlen);
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
            #$nam = substr("$2               ", 0, 15);
            $nam = $2;
        } else {
            die "$ac match\n";
        }


        my $dates = "from $bgn to $end";
        if ($bgn eq $end) {
            #$dates = "on   $bgn               ";
            $dates = "on $bgn";
        }


        if ($org eq "bnbi") {
            $dates = "beginning on $bgn";
        }


        my $str;

        if      ($org eq "craa") {
            #$str .= " *   $nam $dates are Copyright $years Applera Corporation                    and are subject to the GNU General Public License version 2\n";
            $str .= " *   $nam $dates\n";
            $str .= " *     are Copyright $years Applera Corporation, and\n";
            $str .= " *     are subject to the GNU General Public License version 2\n";
            $str .= " *\n";

        } elsif ($org eq "tigr") {
            #$str .= " *   $nam $dates are Copyright $years The Institute for Genomics Research    and are subject to the GNU General Public License version 2\n";
            $str .= " *   $nam $dates\n";
            $str .= " *     are Copyright $years The Institute for Genomics Research, and\n";
            $str .= " *     are subject to the GNU General Public License version 2\n";
            $str .= " *\n";

        } elsif ($org eq "jcvi") {
            #$str .= " *   $nam $dates are Copyright $years J. Craig Venter Institute              and are subject to the GNU General Public License version 2\n";
            $str .= " *   $nam $dates\n";
            $str .= " *     are Copyright $years J. Craig Venter Institute, and\n";
            $str .= " *     are subject to the GNU General Public License version 2\n";
            $str .= " *\n";

        } elsif ($org eq "bnbi") {
            #$str .= " *   $nam $dates are Copyright $years Battelle National Biodefense Institute and are subject to the BSD 3-Clause License\n";
            $str .= " *   $nam $dates\n";
            $str .= " *     are Copyright $years Battelle National Biodefense Institute, and\n";
            $str .= " *     are subject to the BSD 3-Clause License\n";
            $str .= " *\n";

        } elsif ($org eq "nihh") {
            #$str .= " *   $nam $dates are a 'United States Government Work'                       and are released in the public domain\n";
            $str .= " *   $nam $dates\n";
            $str .= " *     are a 'United States Government Work', and\n";
            $str .= " *     are released in the public domain\n";
            $str .= " *\n";

        } elsif ($org eq "none") {
            #$str .= " *   $nam $dates are a 'United States Government Work'                       and are released in the public domain\n";
            $str .= " *   $nam $dates\n";
            $str .= " *     are Copyright $years $nam, and\n";
            $str .= " *     are subject to the GNU General Public License version 2\n";
            $str .= " *\n";

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
    open(L, "> logs");

    for (my $x=7057; $x > 0; $x--) {
        print STDERR "$x\r";

        open(I, "svn log -v -r $x |");
        while (<I>) {
            print L $_;
        }
        close(I);
    }

    close(L);
}



my %authors;
my %copyrights;
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

    #  Deleted a directory, didn't get tracked using the 'from' notation.  Evil.
    $filemap{"src/AS_OVM/overlapInCore-Build_Hash_Index.C"} = "src/overlapInCore/overlapInCore-Build_Hash_Index.C";
    $filemap{"src/AS_OVM/overlapInCore-Extend_Alignment.C"} = "src/overlapInCore/overlapInCore-Extend_Alignment.C";
    $filemap{"src/AS_OVM/overlapInCore-Find_Overlaps.C"} = "src/overlapInCore/overlapInCore-Find_Overlaps.C";
    $filemap{"src/AS_OVM/overlapInCore-Output.C"} = "src/overlapInCore/overlapInCore-Output.C";
    $filemap{"src/AS_OVM/overlapInCore-Process_Overlaps.C"} = "src/overlapInCore/overlapInCore-Process_Overlaps.C";
    $filemap{"src/AS_OVM/overlapInCore-Process_String_Overlaps.C"} = "src/overlapInCore/overlapInCore-Process_String_Overlaps.C";
    $filemap{"src/AS_OVM/overlapInCore-Read_Next_Frag.C"} = "src/overlapInCore/overlapInCore-Read_Next_Frag.C";
    $filemap{"src/AS_OVM/overlapInCore.C"} = "src/overlapInCore/overlapInCore.C";
    $filemap{"src/AS_OVM/overlapInCore.H"} = "src/overlapInCore/overlapInCore.H";
    $filemap{"src/AS_OVM/overlapInCore.mk"} = "src/overlapInCore/overlapInCore.mk";

    #  Deleted a directory, didn't get tracked using the 'from' notation.  Evil.
    $filemap{"src/AS_ALN/AS_ALN_aligners.H"}      = "src/aligners/AS_ALN_aligners.H";
    $filemap{"src/AS_ALN/AS_ALN_bruteforcedp.C"}  = "src/aligners/AS_ALN_bruteforcedp.C";
    $filemap{"src/AS_ALN/AS_ALN_bruteforcedp.H"}  = "src/aligners/AS_ALN_bruteforcedp.H";
    $filemap{"src/AS_ALN/AS_ALN_dpaligner.C"}     = "src/aligners/AS_ALN_dpaligner.C";
    $filemap{"src/AS_ALN/AS_ALN_forcns.C"}        = "src/aligners/AS_ALN_forcns.C";
    $filemap{"src/AS_ALN/AS_ALN_loverlapper.C"}   = "src/aligners/AS_ALN_loverlapper.C";
    $filemap{"src/AS_ALN/AS_ALN_pieceOlap.C"}     = "src/aligners/AS_ALN_pieceOlap.C";
    $filemap{"src/AS_ALN/AS_ALN_qvaligner.C"}     = "src/aligners/AS_ALN_qvaligner.C";
    $filemap{"src/AS_ALN/CA_ALN_local.C"}         = "src/aligners/CA_ALN_local.C";
    $filemap{"src/AS_ALN/CA_ALN_local.H"}         = "src/aligners/CA_ALN_local.H";
    $filemap{"src/AS_ALN/CA_ALN_overlap.C"}       = "src/aligners/CA_ALN_overlap.C";
    $filemap{"src/AS_ALN/aligners.H"}             = "src/aligners/aligners.H";

    #  Branch ovl into ovm, on 2011-07-29
    $filemap{"src/AS_OVL/AS_OVL_overlap_common.h"}  = "src/AS_OVM/overlapInCore-Build_Hash_Index.C";
    $filemap{"src/AS_OVL/AS_OVL_overlap_common.h"} .= "\0src/AS_OVM/overlapInCore-Extend_Alignment.C";
    $filemap{"src/AS_OVL/AS_OVL_overlap_common.h"} .= "\0src/AS_OVM/overlapInCore-Find_Overlaps.C";
    $filemap{"src/AS_OVL/AS_OVL_overlap_common.h"} .= "\0src/AS_OVM/overlapInCore-Output.C";
    $filemap{"src/AS_OVL/AS_OVL_overlap_common.h"} .= "\0src/AS_OVM/overlapInCore-Process_Overlaps.C";
    $filemap{"src/AS_OVL/AS_OVL_overlap_common.h"} .= "\0src/AS_OVM/overlapInCore-Process_String_Overlaps.C";
    $filemap{"src/AS_OVL/AS_OVL_overlap_common.h"} .= "\0src/AS_OVM/overlapInCore-Read_Next_Frag.C";
    $filemap{"src/AS_OVL/AS_OVL_driver_common.h"}   = "src/AS_OVM/overlapInCore.C";
    $filemap{"src/AS_OVL/AS_OVL_overlap.h"}         = "src/AS_OVM/overlapInCore.H";

    $filemap{"src/AS_CNS/MultiAlignment_CNS.h"}     = "src/AS_CNS/MultiAlignment_CNS.h"; 
    $filemap{"src/AS_CNS/MultiAlignment_CNS.h"}    .= "\0src/AS_CNS/MultiAlignment_CNS_private.h"; 

    $filemap{"kmer/leaff/leaff.C"}  = "kmer/leaff/leaff.C";
    $filemap{"kmer/leaff/leaff.C"} .= "kmer/leaff/blocks.C";
    $filemap{"kmer/leaff/leaff.C"} .= "kmer/leaff/dups.C";
    $filemap{"kmer/leaff/leaff.C"} .= "kmer/leaff/gc.C";
    $filemap{"kmer/leaff/leaff.C"} .= "kmer/leaff/partition.C";
    $filemap{"kmer/leaff/leaff.C"} .= "kmer/leaff/stats.C";

    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"}  = "src/AS_CNS/MultiAlignment_CNS.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/AbacusRefine.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/ApplyAlignment.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/BaseCall.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/GetAlignmentTrace.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/MergeMultiAligns.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/MergeRefine.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/MultiAlignContig.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/MultiAlignUnitig.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/PrintAlignment.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/RefreshMANode.c";
    $filemap{"src/AS_CNS/MultiAlignment_CNS.c"} .= "\0src/AS_CNS/ReplaceEndUnitigInContig.c";

    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   = "/src/AS_PER/AS_PER_gkLibrary.C";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/AS_PER_gkStore.C";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/AS_PER_gkStore_UID.C";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/AS_PER_gkStore_clearRange.C";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/AS_PER_gkStream.C";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/gkClearRange.H";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/gkFragment.H";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/gkLibrary.H";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/gkStore.H";
    $filemap{"/src/AS_PER/AS_PER_gkpStore.c"}   .=    "\0src/AS_PER/gkStream.H";

    $filemap{"src/utgcns/libcns/MultiAlignment_CNS.C"}  = "src/utgcns/libcns/abAbacus-populateTig.C";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS.C"} .= "\0src/utgcns/libcns/abAbacus.C";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS.C"} .= "\0src/utgcns/libcns/abColumn.C";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS.C"} .= "\0src/utgcns/libcns/abMultiAlign.C";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS.C"} .= "\0src/utgcns/libcns/abacus-addRead.C";

    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"}  = "src/utgcns/libcns/abAbacus.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abBase.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abBaseCount.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abBead.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abColumn.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abIDs.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abIterators.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abMultiAlign.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abSequence.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/abVariants.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/multiAlignUnitig.H";
    $filemap{"src/utgcns/libcns/MultiAlignment_CNS_private.H"} .= "\0src/utgcns/libcns/unitigConsensus.H";

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
                    $authors{$f}    .= "$A\n";
                    $copyrights{$f} .= "$C\n";
                    $authcopy{$f}   .= "$C$A\n";
                    $derived{$f}    .= "\0$fn";    #  Whatever file we end up with, it was derived from the one we're reading the log for.
                }
            }
        } else {
            print STDERR "WARN: $_\n";
        }
    }
    close(F);
}



foreach my $file (@ARGV) {
    my @lines;

    $file = $1  if ($file =~ m/^\.\/(.*)$/);  #  Remove leading ./ added by find.

    #my @A  = split '\n', $authors{$file};
    #my @C  = split '\n', $copyrights{$file};
    my @AC = split '\n', $authcopy{$file};

    #print STDERR "$file -- a ", scalar(@A), " c ", scalar(@C), "\n";

    #my @Alist  = splitA(@A);
    #my @Clist  = splitC(@C);
    my @AClist = splitAC(@AC);

    my %DElist;
    my @DElist = split '\0', $derived{$file};

    foreach my $d (@DElist) {
        next if ($d eq "");
        #next if ($d eq $file);
        next if (lc $d eq lc $file);
        $DElist{$d}++;
    }

    undef @DElist;

    if (scalar(keys %DElist) > 0) {
        foreach my $d (keys %DElist) {
            push @DElist, " *   $d\n";
        }

        @DElist = sort @DElist;

        unshift @DElist, " *\n";
        unshift @DElist, " * This file is derived from:\n";

        push    @DElist, " *\n";
    }
    


    push @lines, "\n";
    push @lines, "/******************************************************************************\n"; 
    push @lines, " *\n";
    push @lines, " * This file is part of Harly, a software program that assembles whole-genome\n";
    push @lines, " * sequencing reads into contigs.\n";
    push @lines, " *\n";
    push @lines, " * This software is based on RELEASE_1-3_2004-03-17 of the Celera Assembler as\n";
    push @lines, " * distributed by Applera Corporation under the GNU General Public License,\n";
    push @lines, " * version 2.\n";
    push @lines, " *\n";
    push @lines, @DElist;
    push @lines, " * Modifications by:\n";
    push @lines, " *\n";
    push @lines, @AClist;
    push @lines, " * File 'README.licenses' in the root directory of this distribution contains\n";
    push @lines, " * full conditions and disclaimers for each license.\n";
    push @lines, " */\n";
    push @lines, "\n";



    if (0) {
    push @lines, " * For modifications covered by the General Public License (GPL):                                 *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *     This program is free software; you can redistribute it and/or modify it under              *\n";
    push @lines, " *     the terms of the GNU General Public License as published by the Free Software              *\n";
    push @lines, " *     Foundation; either version 2 of the License, or (at your option) any later                 *\n";
    push @lines, " *     version.                                                                                   *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *     This program is distributed in the hope that it will be useful, but WITHOUT ANY            *\n";
    push @lines, " *     WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A            *\n";
    push @lines, " *     PARTICULAR PURPOSE.  See the GNU General Public License for more details.                  *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *     You should have received (LICENSE.txt) a copy of the GNU General Public License            *\n";
    push @lines, " *     along with this program; if not, write to the Free Software Foundation, Inc., 59           *\n";
    push @lines, " *     Temple Place, Suite 330, Boston, MA 02111-1307 USA                                         *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " * For modifications made by the Battelle National Biodefense Institute:                          *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      This Software was prepared for the Department of Homeland Security (DHS) by the           *\n";
    push @lines, " *      Battelle National Biodefense Institute, LLC (BNBI) as part of contract                    *\n";
    push @lines, " *      HSHQDC-07-C-00020 to manage and operate the National Biodefense Analysis and              *\n";
    push @lines, " *      Countermeasures Center (NBACC), a Federally Funded Research and Development               *\n";
    push @lines, " *      Center.                                                                                   *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      Redistribution and use in source and binary forms, with or without modification,          *\n";
    push @lines, " *      are permitted provided that the following conditions are met:                             *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      * Redistributions of source code must retain the above copyright notice, this             *\n";
    push @lines, " *        list of conditions and the following disclaimer.                                        *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      * Redistributions in binary form must reproduce the above copyright notice, this          *\n";
    push @lines, " *        list of conditions and the following disclaimer in the documentation and/or             *\n";
    push @lines, " *        other materials provided with the distribution.                                         *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      * Neither the name of the Battelle National Biodefense Institute nor the names            *\n";
    push @lines, " *        of its contributors may be used to endorse or promote products derived from             *\n";
    push @lines, " *        this software without specific prior written permission.                                *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND           *\n";
    push @lines, " *      ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED             *\n";
    push @lines, " *      WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                    *\n";
    push @lines, " *      DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR          *\n";
    push @lines, " *      ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *\n";
    push @lines, " *      (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;              *\n";
    push @lines, " *      LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON            *\n";
    push @lines, " *      ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                   *\n";
    push @lines, " *      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS             *\n";
    push @lines, " *      SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                              *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " * For modifications made by the National Human Genomics Research Institute:                      *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *                                  PUBLIC DOMAIN NOTICE                                          *\n";
    push @lines, " *                     National Center for Biotechnology Information                              *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      This software/database is a 'United States Government Work' under the terms of            *\n";
    push @lines, " *      the United States Copyright Act.  It was written as part of the author's                  *\n";
    push @lines, " *      official duties as a United States Government employee and thus cannot be                 *\n";
    push @lines, " *      copyrighted.  This software/database is freely available to the public for                *\n";
    push @lines, " *      use. The National Library of Medicine and the U.S.  Government have not placed            *\n";
    push @lines, " *      any restriction on its use or reproduction.                                               *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      Although all reasonable efforts have been taken to ensure the accuracy and                *\n";
    push @lines, " *      reliability of the software and data, the NLM and the U.S.  Government do not             *\n";
    push @lines, " *      and cannot warrant the performance or results that may be obtained by using               *\n";
    push @lines, " *      this software or data. The NLM and the U.S.  Government disclaim all                      *\n";
    push @lines, " *      warranties, express or implied, including warranties of performance,                      *\n";
    push @lines, " *      merchantability or fitness for any particular purpose.                                    *\n";
    push @lines, " *                                                                                                *\n";
    push @lines, " *      Please cite the author in any work or product based on this material.                     *\n";
    }


    #  Scan logs for dates of changes and authors

    @lines = loadFile($file, @lines);

    open(F, "> $file.TEST") or die "Failed to open '$file.TEST' for writing: $!\n";
    print F @lines;
    close(F);
}
