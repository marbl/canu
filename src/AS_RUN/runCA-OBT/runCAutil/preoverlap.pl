use strict;


sub perfectTrimming {
    my $gkpStore = "$wrk/$asm.gkpStore";
    my $refFasta = getGlobal("perfectTrimming");

    return if (!defined($refFasta));

    setGlobal("doOverlapTrimming", 0);

    die "Can't find gkpStore '$gkpStore'\n"  if (! -d $gkpStore);
    die "Can't find reference '$refFasta'\n" if (! -e $refFasta);

    my $cmd;
    my $bin = getBinDirectory();
    my $kmer;
    {
        my @p = split '/', $bin;
        my $l = scalar(@p);

        $p[$l]   = $p[$l-1];
        $p[$l-1] = $p[$l-2];
        $p[$l-2] = "kmer";

        $kmer = join '/', @p;
    }

    if (! -e "$gkpStore/reads.fasta") {
        runCommand($wrk, "$bin/gatekeeper -dumpfastaseq -clear untrim $gkpStore > $gkpStore/reads.fasta") and die;
    }

    if (! -e "$gkpStore/reads.sim4db") {
        #  #1 cov 25 iden 90 sucka
        #  #2 cov 80 iden 96 is nice
        #  #3 cov 25 iden 94 better than #1, but still sucky
        #  #4 cov 80 iden 94 same as #3
        #  #5 cov 80 iden 95
        #  #6 cov 25 iden 96
        #  #7 cov 25 iden 97
        #  #8 cov 25 iden 98
        $cmd  = "$kmer/snapper2";
        $cmd .= " -queries $gkpStore/reads.fasta";
        $cmd .= " -genomic $refFasta";
        $cmd .= " -minhitlength 0";
        $cmd .= " -minhitcoverage 0";
        $cmd .= " -discardexonlength 40";
        $cmd .= " -minmatchcoverage 25";
        $cmd .= " -minmatchidentity 98";
        $cmd .= " -verbose";
        $cmd .= " -numthreads 1";
        $cmd .= " > $gkpStore/reads.sim4db";

        runCommand($wrk, $cmd) and die;
    }

    if (! -e "$gkpStore/reads.extent") {
        runCommand($wrk, "$kmer/pickBestPolish < $gkpStore/reads.sim4db | $kmer/convertToExtent > $gkpStore/reads.extent") and die;
    }

    if (! -e "$gkpStore/reads.update") {
        my %mapBeg;
        my %mapEnd;

        my %allReads;
        my %allMates;

        #  Read the reads and mates
        #
        open(F, "$bin/gatekeeper -dumpfragments -tabular -clear OBT $gkpStore |");
        $_ = <F>;  #  header line
        while (<F>) {
            my @v = split '\s+', $_;
            $allReads{$v[1]}++;
            $allMates{$v[1]} = $v[3];
        }
        close(F);

        #  Read the mapping
        #
        open(F, "< $gkpStore/reads.extent");
        while (<F>) {
            my @v = split '\s+', $_;

            (undef, $v[0]) = split ',', $v[0];

            if ($v[3] < $v[4]) {
                $mapBeg{$v[0]} = $v[3];
                $mapEnd{$v[0]} = $v[4];
            } else {
                $mapBeg{$v[0]} = $v[4];
                $mapEnd{$v[0]} = $v[3];
            }
        }
        close(F);

        #  Update the gkpStore
        #
        open(F, "> $gkpStore/reads.update");
        foreach my $k (keys %allReads) {
            my $mapLen = $mapEnd{$k} - $mapBeg{$k};

            if ($mapLen < 64) {
                print F "frg iid $k isdeleted t\n";
                if ($allMates{$k} > 0) {
                    print F "frg iid $k mateiid 0\n";
                    print F "frg iid $allMates{$k} mateiid 0\n";
                }
            } else {
                print F "frg iid $k orig   $mapBeg{$k} $mapEnd{$k}\n";
                print F "frg iid $k obtini $mapBeg{$k} $mapEnd{$k}\n";
                print F "frg iid $k obt    $mapBeg{$k} $mapEnd{$k}\n";
            }
        }
        close(F);
    }

    if (! -e "$gkpStore/reads.updated") {
        runCommand($wrk, "$bin/gatekeeper -E $gkpStore/reads.update.errors --edit $gkpStore/reads.update $gkpStore > $gkpStore/reads.update.out 2> $gkpStore/reads.update.err") and die;
        touch "$gkpStore/reads.updated";
    }
}


sub preoverlap {
    my @fragFiles = @_;

    $numFrags = getNumberOfFragsInStore($wrk, $asm);

    #  Return if there are fragments in the store, and die if there
    #  are no fragments and no source files.
    #
    if ($numFrags > 0) {
        goto stopafter;
    }

    caFailure("ERROR: No fragment files specified, and stores not already created.\n")
    	if (scalar(@fragFiles) == 0);

    if ((! -d "$wrk/$asm.gkpStore") ||
        (! -e "$wrk/$asm.gkpStore/frg")) {
        my $bin = getBinDirectory();

        #  Make sure all the inputs are here.  We also shred any supplied ace files.
        #
        my $failedFiles = 0;
        my $gkpInput = "";
        foreach my $frg (@fragFiles) {
            if (! -e $frg) {
                print STDERR "MISSING: $frg\n";
                $failedFiles++;
            }

            if ($frg =~ m/^(.*)\.ace$/) {
                my @fff = split '/', $1;
                my $ace = $frg;
                my $nam = pop @fff;

                $frg = "$wrk/$nam.shred.frg";

                if (! -e "$frg") {
                    unlink "$wrk/NB.contigs";
                    unlink "$wrk/NB.shred";

                    if (runCommand($wrk, "perl $bin/Generate_NonShallow_Contigs.pl -a $ace -f $wrk/NB.contigs") ||
                        runCommand($wrk, "perl $bin/Shred_Contigs.pl -f $wrk/NB.contigs > $wrk/NB.shred") ||
                        runCommand($wrk, "perl $bin/FASTA_to_frg_file.pl -f $wrk/NB.shred -q 3 > $frg")) {
                        unlink "$wrk/NB.contigs";
                        unlink "$wrk/NB.shred";
                        unlink "$frg";
                        caFailure("Shredding '$ace' failed.");
                    }

                    unlink "$wrk/NB.contigs";
                    unlink "$wrk/NB.shred";
                }
            }

            if (($frg =~ m/^(.*)\.sff$/) ||
                ($frg =~ m/^(.*)\.sff.gz$/) ||
                ($frg =~ m/^(.*)\.sff.bz2$/)) {
                my @fff = split '/', $1;
                my $sff = $frg;
                my $nam = pop @fff;
                my $log = "$wrk/$nam.sff.log";

                $frg = "$wrk/$nam.sff.frg";

                if (! -e "$frg") {
                    my $bin = getBinDirectory();

                    if (runCommand($wrk, "$bin/sffToCA -libraryname $nam -linker flx -insertsize 3000 300 -log $log -output $frg $sff")) {
                        unlink "$wrk/$frg.sff.frg";
                        caFailure("sffToCA failed.");
                    }
                }
            }

            $gkpInput .= " $frg";
        }
        caFailure("Files supplied on command line not found.\n") if ($failedFiles);

        my $cmd;
        $cmd  = "$bin/gatekeeper ";
        $cmd .= " -o $wrk/$asm.gkpStore.BUILDING ";
        $cmd .= " -T " if (getGlobal("doOverlapTrimming"));
        $cmd .= " -F " if (getGlobal("gkpFixInsertSizes"));
        $cmd .= " -L " if (getGlobal("sffIsPairedEnd") == 1);
        $cmd .= " -E $wrk/$asm.gkpStore.errorLog ";
        $cmd .= "$gkpInput ";
        $cmd .= "> $wrk/$asm.gkpStore.err 2>&1";
        if (runCommand($wrk, $cmd)) {
            caFailure("gatekeeper failed.  Check your input files!\n");
        }

        rename "$wrk/$asm.gkpStore.BUILDING", "$wrk/$asm.gkpStore";
        rmrf("$asm.gkpStore.err");
    }

    perfectTrimming();

    generateVectorTrim();

    my $vi = getGlobal("vectorIntersect");

    if ((defined($vi)) && (! -e "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.log")) {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/gatekeeper -a -v $vi -o $wrk/$asm.gkpStore ";
        $cmd .= "  > $wrk/$asm.gkpStore/$asm.vectorClearLoaded.log";
        $cmd .= " 2> $wrk/$asm.gkpStore/$asm.vectorClearLoaded.err";

        if (runCommand($wrk, $cmd)) {
            rename "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.log", "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.log.FAILED";
            caFailure("Failed.\n");
        }

        rmrf("$wrk/$asm.gkpStore/$asm.vectorClearLoaded.err");
    }

    $numFrags = getNumberOfFragsInStore($wrk, $asm);

  stopafter:
    stopAfter("initialStoreBuilding");
}

1;
