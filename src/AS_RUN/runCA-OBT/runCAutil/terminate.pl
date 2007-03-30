use strict;

sub summarizeConsensusStatistics ($) {
    my $dir = shift @_;

    if (! -e "$dir/consensus.stats.summary") {
        my $NumColumnsInUnitigs           = 0;
        my $NumGapsInUnitigs              = 0;
        my $NumRunsOfGapsInUnitigReads    = 0;
        my $NumColumnsInContigs           = 0;
        my $NumGapsInContigs              = 0;
        my $NumRunsOfGapsInContigReads    = 0;
        my $NumAAMismatches               = 0;
        my $NumFAMismatches               = 0;
        my $NumVARRecords                 = 0;
        my $NumVARStringsWithFlankingGaps = 0;
        my $NumUnitigRetrySuccess         = 0;

        open(F, "ls $dir/$asm*.err |");
        my @files = <F>;
        chomp @files;
        close(F);

        foreach my $f (@files) {
            open(F, "< $f");
            while (<F>) {
                $NumColumnsInUnitigs += $1           if (m/NumColumnsInUnitigs\s+=\s+(\d+)/);
                $NumGapsInUnitigs += $1              if (m/NumGapsInUnitigs\s+=\s+(\d+)/);
                $NumRunsOfGapsInUnitigReads += $1    if (m/NumRunsOfGapsInUnitigReads\s+=\s+(\d+)/);
                $NumColumnsInContigs += $1           if (m/NumColumnsInContigs\s+=\s+(\d+)/);
                $NumGapsInContigs += $1              if (m/NumGapsInContigs\s+=\s+(\d+)/);
                $NumRunsOfGapsInContigReads += $1    if (m/NumRunsOfGapsInContigReads\s+=\s+(\d+)/);
                $NumAAMismatches += $1               if (m/NumAAMismatches\s+=\s+(\d+)/);
                $NumFAMismatches += $1               if (m/NumFAMismatches\s+=\s+(\d+)/);
                $NumVARRecords += $1                 if (m/NumVARRecords\s+=\s+(\d+)/);
                $NumVARStringsWithFlankingGaps += $1 if (m/NumVARStringsWithFlankingGaps\s+=\s+(\d+)/);
                $NumUnitigRetrySuccess += $1         if (m/NumUnitigRetrySuccess\s+=\s+(\d+)/);
            }
            close(F);
        }

        open(F, "> $dir/consensus.stats.summary");
        print F "NumColumnsInUnitigs=$NumColumnsInUnitigs\n"                     if ($NumColumnsInUnitigs > 0);
        print F "NumGapsInUnitigs=$NumGapsInUnitigs\n"                           if ($NumGapsInUnitigs > 0);
        print F "NumRunsOfGapsInUnitigReads=$NumRunsOfGapsInUnitigReads\n"       if ($NumRunsOfGapsInUnitigReads > 0);
        print F "NumColumnsInContigs=$NumColumnsInContigs\n"                     if ($NumColumnsInContigs > 0);
        print F "NumGapsInContigs=$NumGapsInContigs\n"                           if ($NumGapsInContigs > 0);
        print F "NumRunsOfGapsInContigReads=$NumRunsOfGapsInContigReads\n"       if ($NumRunsOfGapsInContigReads > 0);
        print F "NumAAMismatches=$NumAAMismatches\n"                             if ($NumAAMismatches > 0);
        print F "NumFAMismatches=$NumFAMismatches\n"                             if ($NumFAMismatches > 0);
        print F "NumVARRecords=$NumVARRecords\n"                                 if ($NumVARRecords > 0);
        print F "NumVARStringsWithFlankingGaps=$NumVARStringsWithFlankingGaps\n" if ($NumVARStringsWithFlankingGaps > 0);
        print F "NumUnitigRetrySuccess=$NumUnitigRetrySuccess\n"                 if ($NumUnitigRetrySuccess > 0);
        close(F);
    }
}

sub terminate ($) {
    my $cgwDir = shift @_;
    $cgwDir = "$wrk/7-CGW" if (!defined($cgwDir));

    my $termDir = "$wrk/9-terminator";

    if (! -e "$termDir/terminator.success") {
        system("mkdir $termDir") if (! -e "$termDir");

        if (! -e "$termDir/$asm.asm") {
            my $uidServer = getGlobal("uidServer");
            my $fakeUIDs  = getGlobal("fakeUIDs");

            my $cmd;
            $cmd  = "cat $cgwDir/$asm.cgw ";
            $cmd .= " $wrk/8-consensus/$asm.cns_contigs.*[0-9] ";
            $cmd .= " $cgwDir/$asm.cgw_scaffolds | ";
            $cmd .= "$bin/terminator -s $fakeUIDs " if ($fakeUIDs != 0);
            $cmd .= "$bin/terminator -u "           if ($fakeUIDs == 0);
            $cmd .= " $uidServer "                     if (defined($uidServer));
            $cmd .= " -g $wrk/$asm.gkpStore ";
            $cmd .= " -o $termDir/$asm.asm ";
            $cmd .= " -m $termDir/$asm.map ";
            $cmd .= " > $termDir/terminator.err 2>&1 ";

            if (runCommand("$termDir", $cmd)) {
                rename "$termDir/$asm.asm", "$termDir/$asm.asm.FAILED";
                rename "$termDir/$asm.map", "$termDir/$asm.map.FAILED";
                die "Failed.\n";
            }
        }

        ########################################

        if (! -e "$termDir/$asm.scaffold.fasta") {
            my $cmd;
            $cmd  = "$bin/asmProcessScaffolds_TER -q -d ";
            $cmd .= "-f $termDir/$asm.scaffold.fasta ";
            $cmd .= "< $termDir/$asm.asm";

            if (runCommand("$termDir", $cmd)) {
                rename "$termDir/$asm.scaffold.fasta", "$termDir/$asm.scaffold.fasta.FAILED";
                die "Failed.\n";
            }
        }

        ########################################
        #
        #  Generate singletons
        #
        if (! -e "$termDir/$asm.singleton.fasta") {
            my $lastckp = findLastCheckpoint("$wrk/7-CGW");

            my $cmd;
            $cmd  = "$bin/dumpSingletons ";
            $cmd .= " -f $wrk/$asm.frgStore ";
            $cmd .= " -g $wrk/$asm.gkpStore ";
            $cmd .= " -c $cgwDir/$asm -n $lastckp -S ";
            $cmd .= "> $termDir/$asm.singleton.fasta ";
            $cmd .= "2> $termDir/dumpSingletons.err ";

            if (runCommand("$termDir", $cmd)) {
                print STDERR "Failed.\n";
                rename "$termDir/$asm.singleton.fasta", "$termDir/$asm.singleton.fasta.FAILED";
            }
        }

        ########################################

        my $perl = "perl";
        if (-x "/usr/bin/perl") {
            system "/usr/bin/perl -c $bin/caqc.pl >/dev/null 2>&1";
            $perl = "/usr/bin/perl" if ($? == 0);
        }
        if (-x "/usr/local/bin/perl") {
            system "/usr/local/bin/perl -c $bin/caqc.pl >/dev/null 2>&1";
            $perl = "/usr/local/bin/perl" if ($? == 0);
        }


        ########################################
        #
        #  Generate fragment/unitig/contig/scaffold mappings
        #
        if (getGlobal("createPosMap") > 0) {
            if (! -e "$termDir/$asm.posmap.frgscf") {
                my $cmd;
                $cmd = "$perl $bin/buildFragContigPosMap.pl $asm.posmap < $termDir/$asm.asm";
                if (runCommand("$termDir", $cmd)) {
                    rename "$termDir/$asm.posmap.frgscf", "$termDir/$asm.posmap.frgscf.FAILED";
                    die "buildFragContigMap failed.\n";
                }
            }
            if (! -e "$termDir/$asm.posmap.frgscf.sorted") {
                if (runCommand("$termDir", "sort -k2n -T $termDir $termDir/$asm.posmap.frgscf > $termDir/$asm.posmap.frgscf.sorted")) {
                    rename "$termDir/$asm.posmap.frgscf.histogram", "$termDir/$asm.posmap.frgscf.histogram.FAILED";
                } else {
                    runCommand("$termDir", "$bin/fragmentDepth -min       0 -max    3000 < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram1");
                    runCommand("$termDir", "$bin/fragmentDepth -min    3001 -max   10000 < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram2");
                    runCommand("$termDir", "$bin/fragmentDepth -min   10001 -max 1000000 < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram3");
                    runCommand("$termDir", "$bin/fragmentDepth -min 1000001              < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram4");
                }
            }
        }

        ########################################
        #
        #  Sum up the consensus stats, so we can tack them onto the end of
        #  the qc report.
        #
        summarizeConsensusStatistics("$wrk/5-consensus");
        summarizeConsensusStatistics("$wrk/8-consensus");


        ########################################
        #
        #  Generate statistics.
        #
        if (! -e "$termDir/mateLinkIIDRanges.txt") {
            my $cmd;
            $cmd = "$bin/mateLinkIIDRanges.rb $wrk/$asm.gkpStore $bin > $termDir/mateLinkIIDRanges.txt";
            if (runCommand($termDir, $cmd)) {
                rename "$termDir/mateLinkIIDRanges.txt", "$termDir/mateLinkIIDRanges.txt.FAILED";
            }
        }

        if (! -e "$termDir/$asm.qc") {
            if (runCommand("$termDir", "$perl $bin/caqc.pl $termDir/$asm.asm")) {
                rename "$termDir/$asm.qc", "$termDir/$asm.qc.FAILED";
            }

            open(F, ">> $termDir/$asm.qc") or die;

            if (-e "$wrk/5-consensus/consensus.stats.summary") {
                print F "\n[Unitig Consensus]\n";
                open(G, "<  $wrk/5-consensus/consensus.stats.summary") or die;
                while (<G>) {
                    print F $_;
                }
                close(G);
            }

            if (-e "$wrk/8-consensus/consensus.stats.summary") {
                print F "\n[Contig Consensus]\n";
                open(G, "<  $wrk/8-consensus/consensus.stats.summary") or die;
                while (<G>) {
                    print F $_;
                }
                close(G);
            }

            my @H1;
            my @H2;
            my @H3;
            my @H4;
            my $histMax = 0;
            if (-e "$termDir/$asm.posmap.frgscf.histogram1") {
                open(G, "<  $termDir/$asm.posmap.frgscf.histogram1") or die;
                while (<G>) {
                    my ($v, $s) = split '\s+', $_;
                    $H1[$v] = $s;
                    $histMax = $v if ($histMax < $v);
                }
                close(G);
            }
            if (-e "$termDir/$asm.posmap.frgscf.histogram2") {
                open(G, "<  $termDir/$asm.posmap.frgscf.histogram2") or die;
                while (<G>) {
                    my ($v, $s) = split '\s+', $_;
                    $H2[$v] = $s;
                    $histMax = $v if ($histMax < $v);
                }
                close(G);
            }
            if (-e "$termDir/$asm.posmap.frgscf.histogram3") {
                open(G, "<  $termDir/$asm.posmap.frgscf.histogram3") or die;
                while (<G>) {
                    my ($v, $s) = split '\s+', $_;
                    $H3[$v] = $s;
                    $histMax = $v if ($histMax < $v);
                }
                close(G);
            }
            if (-e "$termDir/$asm.posmap.frgscf.histogram4") {
                open(G, "<  $termDir/$asm.posmap.frgscf.histogram4") or die;
                while (<G>) {
                    my ($v, $s) = split '\s+', $_;
                    $H4[$v] = $s;
                    $histMax = $v if ($histMax < $v);
                }
                close(G);
            }

            
            print F "\n[Read Depth Histogram]\n";
            print F "depth    < 3Kbp     < 10Kbp    < 1Mbp     < inf\n";
            for (my $v=0; $v<=$histMax; $v++) {
                printf(F "%-8d %-10d %-10d %-10d %-10d\n",
                       $v, int($H1[$v]), int($H2[$v]), int($H3[$v]), int($H4[$v]));
            }

            close(F);
        }

        if (getGlobal("createAGP") > 0) {
            if (! -e "$termDir/$asm.agp") {
                if (runCommand($termDir, "$perl $bin/asmToAGP.pl < $termDir/$asm.asm > $termDir/$asm.agp")) {
                    rename "$termDir/$asm.agp", "$termDir/$asm.agp.FAILED";
                }
            }
        }

        if (getGlobal("createACE") > 0) {
            if (! -e "$termDir/test.ace.bz2") {
                if (runCommand($termDir, "$bin/gatekeeper -frg $wrk/$asm.gkpStore 2> $termDir/gatekeeper.err | grep -v 'No source' > $termDir/$asm.frg")) {
                    unlink "$termDir/$asm.frg";
                }
                if (runCommand($termDir, "$perl $bin/ca2ace.pl $termDir/$asm.asm")) {
                    rename "$termDir/test.ace.bz2", "$termDir/test.ace.FAILED.bz2";
                }
            }
        }

        touch("$termDir/terminator.success");
        system("ln $termDir/$asm.asm $wrk/$asm.asm");
        system("ln $termDir/$asm.qc  $wrk/$asm.qc");
    }
}

1;
