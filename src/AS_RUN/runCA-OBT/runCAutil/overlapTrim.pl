use strict;

sub overlapTrim {

    return if (getGlobal("doOverlapTrimming") == 0);

    goto alldone if (-e "$wrk/0-overlaptrim/overlaptrim.success");

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");

    #  We use this in a couple of places, so just make it more global.
    #
    my $im = getGlobal("immutableFrags");
    die "immutableFrags '$im' supplied, but not found!\n" if (defined($im) && (! -e $im));

    #  Do an initial overly-permissive quality trimming, intersected
    #  with any known vector trimming.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.initialTrimLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.initialTrimLog.bz2")) {

        my $vi = getGlobal("vectorIntersect");
        my $im = getGlobal("immutableFrags");
        die "vectorIntersect '$vi' supplied, but not found!\n" if (defined($vi) && (! -e $vi));
        die "immutableFrags '$im' supplied, but not found!\n" if (defined($im) && (! -e $im));

        backupFragStore("beforeInitialTrim");

        my $cmd;
        $cmd  = "$bin/initialTrim -update -q 12 ";
        $cmd .= " -vector $vi "    if (defined($vi));
        $cmd .= " -immutable $im " if (defined($im));
        $cmd .= " -log $wrk/0-overlaptrim/$asm.initialTrimLog ";
        $cmd .= " -frg $wrk/$asm.gkpStore ";
        $cmd .= " > $wrk/0-overlaptrim/initialTrim.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            rename "$wrk/0-overlaptrim/$asm.initialTrimLog", "$wrk/0-overlaptrim/$asm.initialTrimLog.failed";
            die "Failed.\n";
        }
    }


    #  Compute overlaps, if we don't have them already

    if (! -e "$wrk/$asm.obtStore") {

        #  Run meryl on the now quality-trimmed frags.  This hopefully
        #  will get around some sequencing centers' habit of having N's in
        #  the low quality region, and is generally just a good idea.
        #
        meryl();

        createOverlapJobs("trim");
        checkOverlap("trim");

        #  Sort the overlaps -- this also duplicates each overlap so that
        #  all overlaps for a fragment A are localized.

        if (runCommand("$wrk/0-overlaptrim",
                       "find $wrk/0-overlaptrim-overlap -follow -name \\*ovb -print > $wrk/0-overlaptrim/all-overlaps-trim.ovllist")) {
            die "Failed to generate a list of all the overlap files.\n";
        }

        my $cmd;
        $cmd  = "$bin/overlapStore ";
        $cmd .= " -c $wrk/$asm.obtStore ";
        $cmd .= " -M " . getGlobal('ovlSortMemory') . " ";
        $cmd .= " -m $numFrags ";
        $cmd .= " -L $wrk/0-overlaptrim/all-overlaps-trim.ovllist";
        $cmd .= " > $wrk/0-overlaptrim/$asm.ovl.sorted.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.ovl.sorted";
            die "Failed to sort.\n";
        }
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated.bz2")) {

        my $cmd;
        $cmd  = "$bin/consolidate ";
        $cmd .= " -ovs $wrk/$asm.obtStore";
        $cmd .= " > $wrk/0-overlaptrim/$asm.ovl.consolidated";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {

          unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated";
          die "Failed to consolidate.\n";
        }
    }


    #  We need to have all the overlaps squashed already, in particular so
    #  that we can get the mode of the 5'mode.  We could do this all in
    #  core, but that would take lots of space.

    if ((! -e "$wrk/0-overlaptrim/$asm.mergeLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.mergeLog.bz2")) {

        backupFragStore("beforeTrimMerge");

        my $cmd;
        $cmd  = "$bin/merge-trimming ";
        $cmd .= "-immutable $im " if (defined($im));
        $cmd .= "-log $wrk/0-overlaptrim/$asm.mergeLog ";
        $cmd .= "-frg $wrk/$asm.gkpStore ";
        $cmd .= "-ovl $wrk/0-overlaptrim/$asm.ovl.consolidated ";
        $cmd .= "> $wrk/0-overlaptrim/$asm.merge.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.mergeLog";
            unlink "$wrk/0-overlaptrim/$asm.mergeLog.stats";
            die "Failed to merge trimming.\n";
        }
    }


    #  Add "-delete" to remove, instead of fix, chimera and spurs.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.chimera.report") &&
        (! -e "$wrk/0-overlaptrim/$asm.chimera.report.bz2")) {

        backupFragStore("beforeChimera");

        my $cmd;
        $cmd  = "$bin/chimera ";
        $cmd .= " -gkp $wrk/$asm.gkpStore ";
        $cmd .= " -ovs $wrk/$asm.obtStore ";
        $cmd .= " -immutable $im " if (defined($im));
        $cmd .= " -summary $wrk/0-overlaptrim/$asm.chimera.summary ";
        $cmd .= " -report  $wrk/0-overlaptrim/$asm.chimera.report ";
        $cmd .= " > $wrk/0-overlaptrim/$asm.chimera.err 2>&1";
        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            rename "$wrk/0-overlaptrim/$asm.chimera.report", "$wrk/0-overlaptrim/$asm.chimera.report.FAILED";
            die "Failed.\n";
        }
    }

    touch("$wrk/0-overlaptrim/overlaptrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}



sub updateFragmentFiles (@) {
    my @frgFiles = @_;

    return if (!defined(getGlobal("updateFragFiles")));

    system("mkdir $wrk/0-overlaptrim-updatedfrg") if (! -d "$wrk/0-overlaptrim-updatedfrg");

    #  Strategy 1) read the list of frag UID and new clear range (and
    #  "is deleted") from the fragstore.  Build a hash here mapping
    #  UID to new ranges.  Stream each fragment file, updating as
    #  needed.

    my %updatedClearRanges;
    my $acc;
    my $deleted;
    open(F, "$bin/dumpFragStore $wrk/$asm.gkpStore | egrep 'Deleted|accID|Orig' |");
    while (<F>) {
        if (m/Deleted:\s*(\d+)\s/) {
            $deleted = $1;
        } elsif (m/accID:(\d+)\s/) {
            $acc = $1;
        } elsif (m/Ovl\((\d+,\d+),\d+\)\s/) {
            if ($deleted ne "1") {
                $updatedClearRanges{$acc} = "$1";
            }
        } else {
            die "Misformed line '$_' updating clear ranges.\n";
        }
    }
    close(F);

    foreach my $f (@frgFiles) {
        my @fn = split '/', $f;
        my $fn = pop @fn;

        if      ($fn =~ m/(.*).gz/) {
            open(F, "gzip -dc $f |");
            $fn = $1;
        } elsif ($fn =~ m/(.*).bz2/) {
            open(F, "bzip2 -dc $f |");
            $fn = $1;
        } else {
            open(F, "< $f");
        }

        my $inFRG;
        my $frg;
        my $inLKG;
        my $lkg;
        my $fg1;
        my $fg2;

        $fn = "$wrk/0-overlaptrim-updatedfrg/$fn";
        open(FN, "> $fn")         or die "Failed to open '$fn'\n";
        open(FD, "> $fn.deleted") or die "Failed to open '$fn.deleted'\n";
        while (<F>) {
            if      (m/^\{FRG$/) {
                $inFRG = 1;
                $frg = $_;
            } elsif ($inFRG) {
                $acc = $1 if (m/acc:(\d+)$/);

                if (m/^clr:\d+,\d+$/) {
                    my $cl = $updatedClearRanges{$acc};
                    if (defined($cl)) {
                        $_ = "clr:$cl\n";
                    }
                }
                $frg .= $_;

                if (m/^\}$/) {
                    if (defined($updatedClearRanges{$acc})) {
                        print FN $frg;
                    } else {
                        print FD $frg;
                    }
                    undef $acc;
                    undef $frg;
                    undef $inFRG;
                }
            } elsif (m/^\{LKG$/) {
                $inLKG = 1;
                $lkg = $_;
            } elsif ($inLKG) {
                $fg1 = $1 if (m/^fg1:(\d+)/);
                $fg2 = $1 if (m/^fg2:(\d+)/);
                $lkg .= $_;

                if (m/^\}$/) {

                    #print STDERR "fg1 -- $fg1 -- $updatedClearRanges{$fg1}\n";
                    #print STDERR "fg2 -- $fg2 -- $updatedClearRanges{$fg2}\n";

                    if (defined($updatedClearRanges{$fg1}) &&
                        defined($updatedClearRanges{$fg2})) {
                        #print STDERR $lkg;
                        print FN $lkg;
                    } else {
                        print FD $lkg;
                    }
                    undef $fg1;
                    undef $fg2;
                    undef $lkg;
                    undef $inLKG;
                }
            } else {
                print FN $_;
            }
        }
        close(FD);
        close(FN);
        close(F);
    }

    print STDERR "Fragment files updated.  Bye.\n";
    exit(0);
}


1;
