use strict;

sub terminate ($) {
    my $cgwDir = shift @_;
    $cgwDir = "$wrk/7-CGW" if (!defined($cgwDir));

    if (! -e "$wrk/$asm.asm") {
        my $uidServer = getGlobal("uidServer");
        my $fakeUIDs  = getGlobal("fakeUIDs");

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "cat $cgwDir/$asm.cgw ";
        $cmd .= " $wrk/8-consensus/$asm.cns_contigs.*[0-9] ";
        $cmd .= " $cgwDir/$asm.cgw_scaffolds | ";
        $cmd .= "$bin/terminator -P -s $fakeUIDs " if ($fakeUIDs != 0);
        $cmd .= "$bin/terminator -P -u "           if ($fakeUIDs == 0);
        $cmd .= " $uidServer "                     if (defined($uidServer));
        $cmd .= " -f $wrk/$asm.frgStore ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -o $wrk/$asm.asm ";
        $cmd .= " -m $wrk/$asm.map ";
        $cmd .= " > $wrk/terminator.err 2>&1 ";

        if (runCommand($cmd)) {
            rename "$wrk/$asm.asm", "$wrk/$asm.asm.FAILED";
            rename "$wrk/$asm.map", "$wrk/$asm.map.FAILED";
            die "Failed.\n";
        }
    }

    ########################################

    if (! -e "$wrk/$asm.scaffold.fasta") {
        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$bin/asmProcessScaffolds_TER -q -d -f $wrk/$asm.scaffold.fasta < $wrk/$asm.asm";
        if (runCommand($cmd)) {
            rename "$wrk/$asm.scaffold.fasta", "$wrk/$asm.scaffold.fasta.FAILED";
            die "Failed.\n";
        }
    }

    ########################################
    #
    #  Generate singletons
    #
    if (! -e "$wrk/$asm.singleton.fasta") {
        my $lastckp = findLastCheckpoint("$wrk/7-CGW");

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$bin/dumpSingletons ";
        $cmd .= " -f $wrk/$asm.frgStore ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -c $cgwDir/$asm -n $lastckp -U ";
        $cmd .= "> $wrk/$asm.singleton.fasta ";
        $cmd .= "2> $wrk/dumpSingletons.err ";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/$asm.singleton.fasta", "$wrk/$asm.singleton.fasta.FAILED";
        }
    }

    ########################################

    my $perl = "perl";
    $perl = "/usr/local/bin/perl" if (-e "/usr/local/bin/perl");
    $perl = "/usr/bin/perl"       if (-e "/usr/bin/perl");

    ########################################
    #
    #  Generate fragment/unitig/contig/scaffold mappings
    #
    if (! -e "$wrk/$asm.posmap.frgscf") {
        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$perl $bin/buildFragContigPosMap.pl $asm.posmap < $wrk/$asm.asm";
        if (runCommand($cmd)) {
            rename "$wrk/$asm.posmap.frgscf", "$wrk/$asm.posmap.frgscf.FAILED";
            die "buildFragContigMap failed.\n";
        }
    }

    ########################################
    #
    #  Sum up the consensus stats, so we can tack them onto the end of
    #  the qc report.
    #
    if (! -e "$wrk/8-consensus/consensus.stats.summary") {
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

        open(F, "ls $wrk/8-consensus/$asm.cns_contigs.*.err |");
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
            }
            close(F);
        }

        open(F, "> $wrk/8-consensus/consensus.stats.summary");
        print F "[Consensus]\n";
        print F "NumColumnsInUnitigs=$NumColumnsInUnitigs\n";
        print F "NumGapsInUnitigs=$NumGapsInUnitigs\n";
        print F "NumRunsOfGapsInUnitigReads=$NumRunsOfGapsInUnitigReads\n";
        print F "NumColumnsInContigs=$NumColumnsInContigs\n";
        print F "NumGapsInContigs=$NumGapsInContigs\n";
        print F "NumRunsOfGapsInContigReads=$NumRunsOfGapsInContigReads\n";
        print F "NumAAMismatches=$NumAAMismatches\n";
        print F "NumFAMismatches=$NumFAMismatches\n";
        print F "NumVARRecords=$NumVARRecords\n";
        print F "NumVARStringsWithFlankingGaps=$NumVARStringsWithFlankingGaps\n";
        close(F);
    }


    ########################################
    #
    #  Generate statistics.
    #
    if (! -e "$wrk/$asm.qc") {

        #  Some amazingly ugly magic to get the perl libs.
        #
        if (!defined($ENV{'PERL5LIB'}) && !defined($ENV{'PERLLIB'})) {
          if (-d "/home/smurphy/preassembly/test/TIGR/scripts") {
            $ENV{'PERL5LIB'} = "/home/smurphy/preassembly/test/TIGR/scripts";
          }
        }

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$perl $bin/caqc.pl $wrk/$asm.asm";
        if (runCommand($cmd)) {
            rename "$wrk/$asm.qc", "$wrk/$asm.qc.FAILED";
            die "Failed.\n";
        }

        open(F, ">> $wrk/$asm.qc") or die;
        print F "\n";

        open(G, "<  $wrk/8-consensus/consensus.stats.summary") or die;
        while (<G>) {
            print F $_;
        }
        close(G);
        close(F);
    }
}

1;
