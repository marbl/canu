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

    if (! -e "$wrk/9-terminator/terminator.success") {
        system("mkdir $wrk/9-terminator") if (! -e "$wrk/9-terminator");

        if (! -e "$wrk/9-terminator/$asm.asm") {
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
            $cmd .= " -o $wrk/9-terminator/$asm.asm ";
            $cmd .= " -m $wrk/9-terminator/$asm.map ";
            $cmd .= " > $wrk/9-terminator/terminator.err 2>&1 ";

            if (runCommand("$wrk/9-terminator", $cmd)) {
                rename "$wrk/9-terminator/$asm.asm", "$wrk/9-terminator/$asm.asm.FAILED";
                rename "$wrk/9-terminator/$asm.map", "$wrk/9-terminator/$asm.map.FAILED";
                die "Failed.\n";
            }
        }

        ########################################

        if (! -e "$wrk/9-terminator/$asm.scaffold.fasta") {
            my $cmd;
            $cmd  = "$bin/asmProcessScaffolds_TER -q -d ";
            $cmd .= "-f $wrk/9-terminator/$asm.scaffold.fasta ";
            $cmd .= "< $wrk/9-terminator/$asm.asm";

            if (runCommand("$wrk/9-terminator", $cmd)) {
                rename "$wrk/9-terminator/$asm.scaffold.fasta", "$wrk/9-terminator/$asm.scaffold.fasta.FAILED";
                die "Failed.\n";
            }
        }

        ########################################
        #
        #  Generate singletons
        #
        if (! -e "$wrk/9-terminator/$asm.singleton.fasta") {
            my $lastckp = findLastCheckpoint("$wrk/7-CGW");

            my $cmd;
            $cmd  = "$bin/dumpSingletons ";
            $cmd .= " -f $wrk/$asm.frgStore ";
            $cmd .= " -g $wrk/$asm.gkpStore ";
            $cmd .= " -c $cgwDir/$asm -n $lastckp -S ";
            $cmd .= "> $wrk/9-terminator/$asm.singleton.fasta ";
            $cmd .= "2> $wrk/9-terminator/dumpSingletons.err ";

            if (runCommand("$wrk/9-terminator", $cmd)) {
                print STDERR "Failed.\n";
                rename "$wrk/9-terminator/$asm.singleton.fasta", "$wrk/9-terminator/$asm.singleton.fasta.FAILED";
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
        if (! -e "$wrk/9-terminator/$asm.posmap.frgscf") {
            my $cmd;
            $cmd = "$perl $bin/buildFragContigPosMap.pl $asm.posmap < $wrk/9-terminator/$asm.asm";
            if (runCommand("$wrk/9-terminator", $cmd)) {
                rename "$wrk/9-terminator/$asm.posmap.frgscf", "$wrk/9-terminator/$asm.posmap.frgscf.FAILED";
                die "buildFragContigMap failed.\n";
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
        my $termDir = "$wrk/9-terminator";
        if (! -e "$termDir/$asm.qc") {

            my $cmd;
            $cmd = "$bin/mateLinkIIDRanges.rb $wrk/$asm.gkpStore $bin > mateLinkIIDRanges.txt";
            if (runCommand($termDir, $cmd)) {
                warn "mateLinkIIDRanges failed mate ranges won't be available.\n";
            }

            $cmd = "$bin/gatekeeper -frg $wrk/$asm.gkpStore > $termDir/$asm.frg 2> $termDir/gatekeeper.err";
            if (runCommand($termDir, $cmd)) {
                warn "gatekeeper didn't dump fragments.\n";
                unlink "$termDir/$asm.frg";
            }
            if (runCommand($termDir, "$perl $bin/ca2ace.pl $asm.asm")) {
                warn "ca2ace failed no ace file created.\n";
            }
            
            if (runCommand($termDir, "$perl $bin/asmToAGP.pl < $asm.asm > $asm.agp")) {
                    warn "asmToAGP.pl failed no agp file created.\n";
            }

            $cmd = "$perl $bin/caqc.pl $wrk/9-terminator/$asm.asm";
            if (runCommand("$wrk/9-terminator", $cmd)) {
                rename "$wrk/9-terminator/$asm.qc", "$wrk/9-terminator/$asm.qc.FAILED";
                die "Failed.\n";
            }

            open(F, ">> $wrk/9-terminator/$asm.qc") or die;
            print F "\n[Unitig Consensus]\n";
            open(G, "<  $wrk/5-consensus/consensus.stats.summary") or die;
            while (<G>) {
                print F $_;
            }
            close(G);

            print F "\n[Contig Consensus]\n";
            open(G, "<  $wrk/8-consensus/consensus.stats.summary") or die;
            while (<G>) {
                print F $_;
            }
            close(G);
            close(F);
        }

        touch("$wrk/9-terminator/terminator.success");
        system("ln $wrk/9-terminator/$asm.asm $wrk/$asm.asm");
        system("ln $wrk/9-terminator/$asm.qc  $wrk/$asm.qc");
    }
}

1;
