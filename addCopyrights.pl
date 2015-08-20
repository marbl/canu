#!/usr/local/bin/perl

use strict;

my @dateStrings = ( "???", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" );




sub getName ($) {
    my $a = $_[0];
    my $A;

    if      ($a eq "adelcher")        {  $A = "Art Delcher";
    } elsif ($a eq "ahalpern")        {  $A = "Aaron Halpern";
    } elsif ($a eq "andreyto")        {
    } elsif ($a eq "florea")          {  $A = "Liliana Florea";
    } elsif ($a eq "cmobarry")        {  $A = "Clark Mobarry";
    } elsif ($a eq "walenz")          {  $A = "Brian P. Walenz";
    } elsif ($a eq "bri")             {  $A = "Brian P. Walenz";
    } elsif ($a eq "brianwalenz")     {  $A = "Brian P. Walenz";
    } elsif ($a eq "catmandew")       {  $A = "Ian Dew";
    } elsif ($a eq "eliv")            {  $A = "Eli Venter";
    } elsif ($a eq "gdenisov")        {  $A = "Gennady Denisov";
    } elsif ($a eq "gesims")          {  $A = "Gregory Sims";
    } elsif ($a eq "granger_sutton")  {  $A = "Granger Sutton";
    } elsif ($a eq "jason_miller")    {  $A = "Jason Miller";
    } elsif ($a eq "jasonmiller9704") {  $A = "Jason Miller";
    } elsif ($a eq "kli1000")         {  $A = "Kelvin Li";
    } elsif ($a eq "mcschatz")        {  $A = "Michael Schatz";
    } elsif ($a eq "mhayton")         {
    } elsif ($a eq "mkotelbajcvi")    {
    } elsif ($a eq "moweis")          {
    } elsif ($a eq "edwardnj")        {  #  kmer build system
    } elsif ($a eq "root")            {  #  Really?  possibly me on os-x
    } elsif ($a eq "halldobv")        {
    } elsif ($a eq "fasulodp")        {  #  kmer build system
    } elsif ($a eq "rbolanos")        {  #  kmer build system
    } elsif ($a eq "ripper")          {  #  kmer build system
    } elsif ($a eq "skoren")          {  $A = "Sergey Koren";
    } elsif ($a eq "vrainish")        {
    } elsif ($a eq "walenzb")         {  $A = "Brian P. Walenz";
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

    if      (($a eq "catmandew") && ($y eq "2004")) {         #  This is the initial commit.

    } elsif  ($a eq "adelcher")        {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "ahalpern")        {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "florea")          {  $C = "none$y$m$d";  #  Earliest 2010-07-08, Latest 2011-11-16
    } elsif  ($a eq "cmobarry")        {  $C = "craa$y$m$d";
    } elsif  ($a eq "catmandew")       {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "eliv")            {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "gdenisov")        {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "gesims")          {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "granger_sutton")  {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "jason_miller")    {  $C = "tigr$y$m$d";
    } elsif  ($a eq "jasonmiller9704") {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "kli1000")         {  $C = "jcvi$y$m$d";
    } elsif  ($a eq "mcschatz")        {  $C = "tigr$y$m$d";

    } elsif (($a eq "skoren") && ($y <  "2011")) {  $C = "jcvi$y$m$d";
    } elsif (($a eq "skoren") && ($y >= "2011")) {  $C = "bnbi$y$m$d";

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

        if ($org eq "bnbi") {
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

    for (my $x=$rev; $x > 0; $x--) {
        print STDERR "$x\r";

        open(I, "svn log -v -r $x |");
        while (<I>) {
            print L $_;
        }
        close(I);
    }

    close(L);
}

#  Build a mapping from original file name -> new file name.  Because I didn't use svn copy to
#  branch files, and because some files exploded into smaller pieces.

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

    #  Deleted a directory, didn't get tracked using the 'from' notation.  Evil.
    $filemap{"src/pipelines/ca3g/Consensus.pm"}              = "src/pipelines/canu/Consensus.pm";
    $filemap{"src/pipelines/ca3g/CorrectReads.pm"}           = "src/pipelines/canu/CorrectReads.pm";
    $filemap{"src/pipelines/ca3g/Defaults.pm"}               = "src/pipelines/canu/Defaults.pm";
    $filemap{"src/pipelines/ca3g/Execution.pm"}              = "src/pipelines/canu/Execution.pm";
    $filemap{"src/pipelines/ca3g/Gatekeeper.pm"}             = "src/pipelines/canu/Gatekeeper.pm";
    $filemap{"src/pipelines/ca3g/Meryl.pm"}                  = "src/pipelines/canu/Meryl.pm";
    $filemap{"src/pipelines/ca3g/Output.pm"}                 = "src/pipelines/canu/Output.pm";
    $filemap{"src/pipelines/ca3g/OverlapBasedTrimming.pm"}   = "src/pipelines/canu/OverlapBasedTrimming.pm";
    $filemap{"src/pipelines/ca3g/OverlapErrorAdjustment.pm"} = "src/pipelines/canu/OverlapErrorAdjustment.pm";
    $filemap{"src/pipelines/ca3g/OverlapInCore.pm"}          = "src/pipelines/canu/OverlapInCore.pm";
    $filemap{"src/pipelines/ca3g/OverlapMhap.pm"}            = "src/pipelines/canu/OverlapMhap.pm";
    $filemap{"src/pipelines/ca3g/OverlapStore.pm"}           = "src/pipelines/canu/OverlapStore.pm";
    $filemap{"src/pipelines/ca3g/Unitig.pm"}                 = "src/pipelines/canu/Unitig.pm";

    #  Just a tiny bit of the executive survived, but it's basically in every file.
    #  Execution.pm has a significant chunk from ESTmapper/scheduler.pm (and then in runCA.pl)
    $filemap{"src/AS_RUN/runCA.pl"}          = "src/pipelines/canu/Execution.pl";
    $filemap{"kmer/ESTmapper/scheduler.pm"} .= "\0src/pipelines/canu/Execution.pm";
    $filemap{"kmer/scripts/libBri.pm"}      .= "\0src/pipelines/canu/Execution.pm";

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
                    $authcopy{$f}   .= "$C$A\n";
                    $derived{$f}    .= "$fn\n";    #  Whatever file we end up with, it was derived from the one we're reading the log for.
                }
            }
        } else {
            print STDERR "WARN: $_\n";
        }
    }
    close(F);
}





#  Process each file.

open(FIN, "find kmer src -type f -print |") or die "Failed to launch 'find'\n";
open(OUT, "> addCopyrights.dat") or die "Failed to open 'addCopyrights.dat' for writing: $!\n";

while (<FIN>) {
    chomp;

    my $file = $_;

    $file = $1  if ($_ =~ m/^\.\/(.*)$/);  #  Remove leading ./ added by find.

    my @lines;


    next if ($file =~ m/\.mk$/);
    next if ($file =~ m/Makefile/);

    next if ($file =~ m/\.sh$/);

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
            print OUT "D\t$file\t$d\n";
        }

        @DElist = sort @DElist;

        unshift @DElist, " $cc\n";
        unshift @DElist, " $cc  This file is derived from:\n";

        push    @DElist, " $cc\n";
    }

    if ($file =~ m/\.pl$/) {
        push @lines, "#!perl\n";
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
        s/\s+$//;

        #  If not at the start, add the line.

        if ($start == 0) {
            push @lines, "$_\n";
            next;
        }

        #  Else, we're at the start; if blank or a comment, skip it.  Only C-style comments are skipped.

        if (($_ eq "") || ($_ =~ m/^[\/\s]\*/) || (($_ =~ m/^\s*#/) && ($_ !~ m/^#include/))) {
            next;
        }

        #  Else, add the line, and declare that we're no longer at the start.

        push @lines, "$_\n";

        $start = 0;
    }
    close(F);

    rename "$file", "$file.ORIG";

    open(F, "> $file") or die "Failed to open '$file' for writing: $!\n";
    print F @lines;
    close(F);

    #open(F, "> $file.MODIFIED") or die "Failed to open '$file.MODIFIED' for writing: $!\n";
    #print F @lines;
    #close(F);
}

close(FIN);
