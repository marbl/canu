use strict;

my $terminateFakeUID = 0;


sub terminate {

    my $failedJobs = 0;

    open(CGWIN, "ls $wrk/7-CGW/$asm.cgw_contigs.* |") or die;
    while (<CGWIN>) {
        chomp;

        if (m/cgw_contigs.(\d+)/) {
            my $jobName   = substr("0000000000" . $1, -8);

            if (! -e "$wrk/8-consensus/$jobName.success") {
                print STDERR "$wrk/8-consensus/$jobName failed.\n";
                $failedJobs++;
            }
        } else {
            print STDERR "WARNING: didn't match $_ for cgw_contigs filename!\n";
        }
    }
    close(CGWIN);

    if ($failedJobs) {
        print STDERR "$failedJobs failed.  Good luck.\n";
        exit(1);
    }

    ########################################

    if (! -e "$wrk/$asm.asm") {
        print STDERR "Starting h -- terminator\n";

        my $uidServer = getGlobal("uidServer", undef);

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "cat $wrk/7-CGW/$asm.cgw ";
        $cmd .= "    $wrk/8-consensus/$asm.cns_contigs.*[0-9] ";
        $cmd .= "    $wrk/7-CGW/$asm.cgw_scaffolds | ";
        $cmd .= "$bin/terminator -P -s 100000000 " if ($terminateFakeUID != 0);
        $cmd .= "$bin/terminator -P -u "           if ($terminateFakeUID == 0);
        $cmd .= "    $uidServer "                  if (defined($uidServer));
        $cmd .= "    -f $wrk/$asm.frgStore ";
        $cmd .= "    -g $wrk/$asm.gkpStore ";
        $cmd .= "    -o $wrk/$asm.asm ";
        $cmd .= "    -m $wrk/$asm.map ";
        $cmd .= "> $wrk/terminator.out ";
        $cmd .= "2> $wrk/terminator.err";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/$asm.asm", "$wrk/$asm.asm.FAILED";
            rename "$wrk/$asm.map", "$wrk/$asm.map.FAILED";
            exit(1);
        }
    }

    ########################################

    if (! -e "$wrk/$asm.scaffold.fasta") {
        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$bin/asmProcessScaffolds_TER -q -d -f $wrk/$asm.scaffold.fasta < $wrk/$asm.asm";
        if (runCommand($cmd)) {
            print "Failed.\n";
            rename "$wrk/$asm.scaffold.fasta", "$wrk/$asm.scaffold.fasta.FAILED";
            exit(1);
        }
    }

    ########################################
    #
    #  Generate singletons
    #
    # /bioinfo/work/projects/macaque-v3/wgs/Linux64/bin/dumpSingletons
    # -f v1.frgStore -g v1.gkpStore -c 7-4-CGW/v1 -n 85 -U > v1.singleton.fasta



    ########################################

    my $perl = "perl";
    $perl = "/usr/local/bin/perl" if (-e "/usr/local/bin/perl");
    $perl = "/usr/bin/perl"       if (-e "/usr/bin/perl");

    if (0) {
        if (! -e "$wrk/$asm.scflen") {
            my $cmd;
            $cmd  = "cd $wrk && ";
            $cmd .= "$perl /home/ahalpern/asm_parse.pl $wrk/$asm.frgctg $wrk/$asm.ctglen $wrk/$asm.ctgscf $wrk/$asm.scflen ";
            $cmd .= "$wrk/$asm.asm";
            if (runCommand($cmd)) {
                print "Failed.\n";
                rename "$wrk/$asm.frgctg", "$wrk/$asm.frgctg.FAILED";
                rename "$wrk/$asm.ctglen", "$wrk/$asm.ctglen.FAILED";
                rename "$wrk/$asm.ctgscf", "$wrk/$asm.ctgscf.FAILED";
                rename "$wrk/$asm.scflen", "$wrk/$asm.scflen.FAILED";
                exit(1);
            }
        }

        if (! -e "$wrk/$asm.frgscf") {
            my $cmd;
            $cmd  = "cd $wrk && ";
            $cmd .= "$perl /home/ahalpern/frg_onto_scf_map.pl $wrk/$asm.ctgscf $wrk/$asm.frgctg > $wrk/$asm.frgscf";
            if (runCommand($cmd)) {
                print "Failed.\n";
                rename "$wrk/$asm.frgscf", "$wrk/$asm.frgscf.FAILED";
                exit(1);
            }
        }
    }

    ########################################

    #  Generate statistics.  There be magic here.
    #  It lives in CVS under tools/asm_scripts.
    #  But not in the assembly tree.

    if (! -e "$wrk/$asm.qc") {
        $ENV{'PERL5LIB'} = "/home/smurphy/preassembly/test/TIGR/scripts";

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$perl /jtc/datadelivery/bin/caqc.pl $wrk/$asm.asm";
        if (runCommand($cmd)) {
            print "Failed.\n";
            rename "$wrk/$asm.qc", "$wrk/$asm.qc.FAILED";
            exit(1);
        }
    }
}

1;
