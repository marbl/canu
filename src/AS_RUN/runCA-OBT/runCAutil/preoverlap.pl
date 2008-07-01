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

    system("mkdir $wrk/0-preoverlap") if (! -d "$wrk/0-preoverlap");

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

                my $w = "$wrk/0-preoverlap";

                $frg = pop @fff;
                $frg = "$w/$frg.frg";

                unlink "$w/NB.contigs";
                unlink "$w/NB.shred";
                runCommand($w, "perl $bin/Generate_NonShallow_Contigs.pl -a $ace -f $w/NB.contigs");
                runCommand($w, "perl $bin/Shred_Contigs.pl -f $w/NB.contigs > $w/NB.shred");
                runCommand($w, "perl $bin/FASTA_to_frg_file.pl -f $w/NB.shred -q 3 > $frg");
                unlink "$w/NB.contigs";
                unlink "$w/NB.shred";
            }

            $gkpInput .= " $frg";
        }
        caFailure("Files supplied on command line not found.\n") if ($failedFiles);

        my $cmd;
        $cmd  = "$bin/gatekeeper ";
        $cmd .= " -o $wrk/$asm.gkpStore ";
        $cmd .= " -T " if (getGlobal("doOverlapTrimming"));
        $cmd .= " -F " if (getGlobal("gkpFixInsertSizes"));
        $cmd .= " -L " if (getGlobal("sffIsPairedEnd") == 1);
        $cmd .= " -E $wrk/0-preoverlap/gatekeeper.errors ";
        $cmd .= "$gkpInput ";
        $cmd .= "> $wrk/0-preoverlap/gatekeeper.err 2>&1";
        if (runCommand("$wrk/0-preoverlap", $cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/0-preoverlap/$asm.inp", "$wrk/0-preoverlap/$asm.inp.FAILED";
            rename "$wrk/$asm.gkpStore", "$wrk/$asm.gkpStore.FAILED";
            caFailure("gatekeeper failed.  Check your input files!\n");
        }
    }

    perfectTrimming();

    generateVectorTrim();
    my $vi = getGlobal("vectorIntersect");
    if ((defined($vi)) && (! -e "$wrk/0-preoverlap/$asm.vectorClearLoaded")) {
        my $bin = getBinDirectory();
        if (runCommand("$wrk/0-preoverlap", "$bin/gatekeeper -a -v $vi -o $wrk/$asm.gkpStore > $wrk/0-preoverlap/$asm.vectorClearLoaded.err 2>&1")) {
            caFailure("Failed.\n");
        }
        touch("$wrk/0-preoverlap/$asm.vectorClearLoaded");
    }

    $numFrags = getNumberOfFragsInStore($wrk, $asm);

    #if ( ! -s "$wrk/$asm.frg" ) { # don't overwrite if it's already there
    #    my $bin = getBinDirectory();
    #    if (runCommand($wrk, "$bin/gatekeeper -dumpfrg $wrk/$asm.gkpStore 2> $wrk/gatekeeper.err | grep -v 'No source' > $wrk/$asm.frg")) {
    #        unlink "$wrk/$asm.frg";
    #    }
    #}

  stopafter:
    stopAfter("initialStoreBuilding");
}

1;
