use strict;

sub terminate ($) {
    my $cgwDir = shift @_;
    $cgwDir = "$wrk/7-CGW" if (!defined($cgwDir));

    my $failedJobs       = 0;

    open(CGWIN, "ls $cgwDir/$asm.cgw_contigs.* |") or die;
    while (<CGWIN>) {
        chomp;

        if (m/cgw_contigs.(\d+)/) {
            if (! -e "$wrk/8-consensus/$asm.cns_contigs.$1.success") {
                print STDERR "$wrk/8-consensus/$asm.cns_contigs.$1 failed.\n";
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
        $cmd .= "> $wrk/terminator.err 2>&1 ";

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
    if (! -e "$wrk/$asm.singleton.fasta") {
        my $lastckp = findLastCheckpoint("$wrk/7-CGW");

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$bin/dumpSingletons ";
        $cmd .= " -f $wrk/$asm.frgStore ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -c $cgwDir/$asm -n $lastckp -U > $wrk/$asm.singleton.fasta";
        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/$asm.singleton.fasta", "$wrk/$asm.singleton.fasta.FAILED";
        }
    }

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
        if (!defined($ENV{'PERL5LIB'})) {
          if (-d "/home/smurphy/preassembly/test/TIGR/scripts") {
            $ENV{'PERL5LIB'} = "/home/smurphy/preassembly/test/TIGR/scripts";
          } elsif (defined $ENV{PERLLIB}) {
            ## MCS: For some reason at CBCB PERL5LIB needs to be set from PERLLIB
            $ENV{'PERL5LIB'} = $ENV{PERLLIB};
          }
        }

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$perl $bin/caqc.pl $wrk/$asm.asm";
        if (runCommand($cmd)) {
            print "Failed.\n";
            rename "$wrk/$asm.qc", "$wrk/$asm.qc.FAILED";
            exit(1);
        }
    }
}

1;
