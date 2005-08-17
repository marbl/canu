#!/usr/bin/perl

use strict;
use Config;  #  for @signame

#  Does an assembly.  Optionally does overlap-based trimming.

my $merThresholdT = 100;  #  For trimming
my $merThresholdA =  50;  #  For assembly

my $doTrimming           = 1;
my $doVectorIntersection = undef;
my $doFixChimera         = 1;

my $meryl     = "/bioinfo/assembly/walenz/src/genomics/meryl/meryl";
my $tbin      = "/bioinfo/assembly/walenz/trim/trimming";
my $abin      = "/bioinfo/assembly/walenz/WGA/wgs-assembler/OSF1/bin";
my $prefix    = undef;

my $uname = `uname -m`;
if ($uname =~ m/x86_64/) {
    $meryl     = "/bioinfo/work/projects/walenz/src/genomics/meryl/meryl";
    $tbin      = "/bioinfo/work/projects/walenz/trim/trimming";
    $meryl     = "/project/dros/dros-clear-trim/src/genomics/meryl/meryl";
    $tbin      = "/project/dros/dros-clear-trim/trimming";
    $abin      = "/bioinfo/assembly/walenz/WGA/wgs-assembler/Linux64/bin";
    #$abin      = "/project/dros/dros-clear-trim/wga/Linux64/bin";
    #print STDERR "USING HACKED WGA SOURCE!\n";
}
if ($uname =~ m/i686/) {
    #  Not really supported -- you have to rebuild src/genomics/meryl and copy it
    #  into the asm bin directory.
    #
    $meryl     = "/bioinfo/assembly/walenz/src/genomics/meryl/meryl";
    $meryl     = "/bioinfo/assembly/walenz/WGA/wgs-assembler/Linux/bin/meryl2";
    $tbin      = "/bioinfo/assembly/walenz/trim/trimming";
    $abin      = "/bioinfo/assembly/walenz/WGA/wgs-assembler/Linux/bin";
}

#$meryl = "/home/work/src/genomics/meryl/meryl";
#$tbin  = "/home/work/trim/trimming";
#$abin  = "/home/work/WGA/wgs-assembler/Linux/bin";

print STDERR "WARNING: USING -f TO CONSENSUS!\n";

if (-e "doasm.opts") {
    print STDERR "Options file found, ignoring command line options!\n";

    open(F, "< doasm.opts");
    $prefix               = <F>;  chomp $prefix;
    $doTrimming           = <F>;  chomp $doTrimming;
    $merThresholdT        = <F>;  chomp $merThresholdT;
    $merThresholdA        = <F>;  chomp $merThresholdA;
    $doVectorIntersection = <F>;  chomp $doVectorIntersection;
    $doFixChimera         = <F>;  chomp $doFixChimera;
    close(F);
} else {
    while (scalar(@ARGV)) {
        my $arg = shift @ARGV;
        if      ($arg eq "-p") {
            $prefix = shift @ARGV;
        } elsif ($arg eq "-notrim") {
            $doTrimming = 0;
        } elsif ($arg eq "-trim") {
            $doTrimming = 1;
        } elsif ($arg eq "-vectorintersection") {
            $doVectorIntersection = shift @ARGV;
        } elsif ($arg eq "-fixchimera") {
            $doFixChimera = shift @ARGV;
        } elsif ($arg eq "-mttrim") {
            $merThresholdT = int(shift @ARGV);
        } elsif ($arg eq "-mtasm") {
            $merThresholdA = int(shift @ARGV);
        } else {
            print STDERR "UNKNOWN OPTION $arg\n";
        }
    }
}

if (!defined($prefix)) {
    print STDERR "usage: $0 -p prefix\n";
    exit(1);
}

open(F, "> doasm.opts");
print F "$prefix\n";
print F "$doTrimming\n";
print F "$merThresholdT\n";
print F "$merThresholdA\n";
print F "$doVectorIntersection\n";
print F "$doFixChimera\n";
close(F);


#  Build a new gkpStore, frgStore.
#
if (! -e "$prefix.gkpStore") {
    runCommand("$abin/gatekeeper -e 999999 -X -C -N -Q -P -f $prefix.gkpStore $prefix.frg") and die;
}
if (! -e "$prefix.frgStore") {
    runCommand("$abin/PopulateFragStore -P -c -f -o $prefix.frgStore -V $prefix.ofg $prefix.inp") and die;
}


#  We use the latest meryl, but the CA meryl would also work.  Two
#  sets of mers are generated, one at the merThreshold for use in the
#  assembly overlap phase, and one at twice the merThreshold for use
#  in clear range determination.
#
if (! -e "$prefix.nmers$merThresholdA.fasta") {
    if ((! -e "$prefix.mcidx") || (! -e "$prefix.mcdat")) {
        runCommand("cat $prefix.frg | $meryl -v -B -C -m 22 -s - -o $prefix") and die;
    }
    runCommand("$meryl -v -Dt -n $merThresholdT -s $prefix > $prefix.nmers$merThresholdT.fasta") and die;
    runCommand("$meryl -v -Dt -n $merThresholdA -s $prefix > $prefix.nmers$merThresholdA.fasta") and die;
}

if ($doTrimming) {

    #  Do a leniant quality filter.  Run overlapper with the Granger
    #  option (-G).  We used to fiddle with the sequences to convert
    #  any N into a random base with low quality.

    if (! -e "$prefix.trim.qualityLog") {
        if (runCommand("$tbin/qualityTrim -update -log $prefix.trim.qualityLog -q 12 -frg $prefix.frgStore")) {
            rename "$prefix.trim.quailtyLog", "$prefix.trim.qualityLog.failed";
            die "Failed.\n";
        }
    }

    #  Do the _optional_ vector intersection

    if ($doVectorIntersection) {
        if (! -e "$prefix.trim.vectorIntersectionLog") {
            if (runCommand("$tbin/intersectTrim -update -intersect $doVectorIntersection -log $prefix.trim.vectorIntersectionLog -frg $prefix.frgStore")) {
                rename "$prefix.trim.vectorIntersectionLog", "$prefix.trim.vectorIntersectionLog.failed";
                die "Failed.\n";
            }
        }
    }


    if (! -e "$prefix.trim.ovl") {
        if (runCommand("$abin/overlap -M 4GB -G -t 4 -P -h 1- -r 1- -k $prefix.nmers$merThresholdT.fasta -o $prefix.trim.ovl $prefix.frgStore")) {
            rename "$prefix.trim.ovl", "$prefix.trim.ovl.failed";
            die "Failed.\n";
        }
    }

    #  Sort the overlaps -- this also duplicates each overlap so that
    #  all overlaps for a fragment A are localized.

    if (! -e "$prefix.trim.ovl.sorted") {
      my $maxiid = 0;
      open(F, "$abin/lastfraginstore $prefix.frgStore |") or die "Failed to lastfraginstore.";
      while (<F>) {
        if (m/Last frag in store is iid = (\d+)/) {
          $maxiid = $1;
        }
      }
      close(F);

      if ($maxiid == 0) {
        die "Failed to find the number of frags in the store!\n";
      }

      if (runCommand("$tbin/sort-overlaps -memory 16000 -maxiid $maxiid $prefix.trim.ovl > $prefix.trim.ovl.sorted")) {
        unlink "$prefix.trim.ovl.sorted";
        die "Failed to sort.\n";
      }
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if (! -e "$prefix.trim.ovl.consolidated") {
        if (runCommand("$tbin/consolidate < $prefix.trim.ovl.sorted > $prefix.trim.ovl.consolidated.1")) {
          unlink "$prefix.ovl.trim.consolidated.1";
          die "Failed to sort.\n";
        }
       
        #  Clean up stuff
        #   - add missing fragments to $prefix.trim.ovl.consolidated
        #
        open(F, "< $prefix.trim.ovl.consolidated.1");
        open(G, "> $prefix.trim.ovl.consolidated");
        my $inId = 0;
        my $otId = 0;
        while (<F>) {
            ($inId) = split '\s+', $_;
            $otId++;
            while ($otId < $inId) {
                #print STDERR "$otId has no overlaps (but $inId does).\n";
                print G "$otId  0 0 0 0 0  0 0 0 0 0  0\n";
                $otId++;
            }
            print G $_;
            $otId = $inId;
        }
        close(G);
        close(F);
    }


    #  We need to have all the overlaps squashed already, in particular so
    #  that we can get the mode of the 5'mode.  We could do this all in
    #  core, but that would take lots of space.

    #  This is for restarting -- I always seem to remove the *.ofg, and
    #  forget to rename the original back.
    if ((! -e "$prefix.ofg") && (-e "$prefix.ofg.orig")) {
        rename "$prefix.ofg.orig", "$prefix.ofg";
    }

    if (! -e "$prefix.ofg.orig") {
        runCommand("$tbin/merge-trimming -log $prefix.trim.mergeLog -frg $prefix.frgStore -ovl $prefix.trim.ovl.consolidated") and die;
        runCommand("$abin/dumpFragStoreAsOFG $prefix.frgStore > $prefix.2.ofg") and die;
        rename "$prefix.ofg", "$prefix.ofg.orig";
        rename "$prefix.2.ofg", "$prefix.ofg";
    }


    #  Be nice, and generate a report for Granger
    #
    if (! -e "$prefix.trim.report") {
        open(A, "< $prefix.trim.qualityLog") or die "Failed to open $prefix.trim.qualityLog\n";
        open(B, "< $prefix.trim.mergeLog") or die "Failed to open $prefix.trim.mergeLog\n";
        open(C, "< $prefix.trim.ovl.consolidated") or die "Failed to open $prefix.trim.ovl.consolidated\n";
        open(F, "> $prefix.trim.report") or die "Failed to open $prefix.trim.report\n";

        while (!eof(A) || !eof(B) || !eof(C)) {
            my $a = <A>; chomp $a;
            my $b = <B>; chomp $b;
            my $c = <C>; chomp $c;

            my @av = split '\s+', $a;
            my @bv = split '\s+', $b;
            my @cv = split '\s+', $c;

            if (($av[0] != $bv[0]) || ($bv[0] != $cv[0]) || ($av[0] != $cv[0])) {
                print STDERR "ERROR: ID MISMATCH!\n";
                print STDERR "A: $a\nB: $b\nC: $c\n";
            }

            printf(F "%6d : TI: %4d %4d Q1: %4d %4d Q2: %4d %4d TF: %4d %4d : %s\n",
                   $av[0],
                   $av[1], $av[2],  #  TI
                   $av[4], $av[5],  #  Q1
                   $bv[1], $bv[2],  #  Q2
                   $bv[3], $bv[4],  #  TF
                   $c);
        }

        close(C);
        close(B);
        close(A);
        close(F);
    }

    if (! -e "$prefix.trim.chimera.report") {
        my $delete = "-delete";
        if ($doFixChimera) {
          $delete = "";
        }
        runCommand("$tbin/chimera $delete -frg $prefix.frgStore < $prefix.trim.ovl.sorted > $prefix.trim.chimera.report") and die;
    }

}


if (! -e "$prefix.ovl") {
    if (runCommand("$abin/overlap -M 4GB -t 4 -P -h 1- -r 1- -k $prefix.nmers$merThresholdA.fasta -o $prefix.ovl $prefix.frgStore")) {
        rename "$prefix.ovl", "$prefix.ovl.FAILED";
        die "Failed.\n";
    }
}


if (! -e "$prefix.ovlStore") {
    open(F, "> $prefix.ovllist");
    print F "$prefix.ovl\n";
    close(F);
    if (runCommand("$abin/grow-olap-store -M 4096 -cfS -o $prefix.ovlStore -L $prefix.ovllist")) {
        rename "$prefix.ovlStore", "$prefix.ovlStore.FAILED";
        die "Failed.\n";
    }
    unlink "$prefix.ovllist";
}


open(F, "$abin/lastfraginstore $prefix.frgStore |");
$_ = <F>;
chomp $_;
close(F);
my $fend = 0;
if (m/^Last frag in store is iid = (\d+)$/) {
    $fend = $1;
} else {
    print STDERR "Failed to determine the last frag in store!  Got '$_'\n";
    exit(1);
}

if (! -e "$prefix.1.corr") {
    if ($fend > 1000000) {
	my $fb = 1;
        my $fe = int($fend / 3);
	my $fs = $fe;

        if (runCommand("$abin/correct-frags -t 4 -k 9 -S $prefix.ovlStore -x 1 -o $prefix.1.corr $prefix.frgStore $fb $fe")) {
            rename "$prefix.1.corr", "$prefix.1.corr.FAILED";
            die "Failed.\n";
        }
        $fb  = $fe+1;
        $fe += $fs;
        if (runCommand("$abin/correct-frags -t 4 -k 9 -S $prefix.ovlStore -x 1 -o $prefix.2.corr $prefix.frgStore $fb $fe")) {
            rename "$prefix.2.corr", "$prefix.2.corr.FAILED";
            die "Failed.\n";
        }
        $fb  = $fe+1;
        $fe  = $fend;
        if (runCommand("$abin/correct-frags -t 4 -k 9 -S $prefix.ovlStore -x 1 -o $prefix.3.corr $prefix.frgStore $fb $fe")) {
            rename "$prefix.3.corr", "$prefix.3.corr.FAILED";
            die "Failed.\n";
        }
    } else {
        if (runCommand("$abin/correct-frags -t 4 -k 9 -S $prefix.ovlStore -x 1 -o $prefix.1.corr $prefix.frgStore 1 $fend")) {
            rename "$prefix.1.corr", "$prefix.1.corr.FAILED";
            die "Failed.\n";
        }
    }
}


if (! -e "$prefix.erate") {
    system("ls -1 $prefix.*.corr > corrlist");
    if (runCommand("$abin/cat-corrects -o $prefix.corr -L corrlist")) {
        rename "$prefix.corr", "$prefix.corr.FAILED";
        die "Failed\n";
    }
    unlink "corrlist";

    if (runCommand("$abin/correct-olaps -S $prefix.ovlStore -e $prefix.erate $prefix.frgStore $prefix.corr 1 $fend")) {
        rename "$prefix.erate", "$prefix.erate.FAILED";
        die "Failed.\n";
    }
}


if (! -e "$prefix.erate.updated") {
    if (runCommand("$abin/update-erates $prefix.ovlStore $prefix.erate")) {
        die "Failed.\n";
    }
    open(F, "> $prefix.erate.updated");
    close(F);
}


if (! -e "$prefix.cgb") {
    #  -U 1 -- BUBBLESMOOTHING
    #  -m PREEDGES     -- preallocate memory
    #  -e OVERLAPERROR -- errors, 1.5%, overlaps with more than this are discarded
    #
    open(F, "> ofglist");
    print F "$prefix.ofg\n";
    close(F);
    if (runCommand("$abin/unitigger -c -P -A 1 -U 1 -e 15 -n $fend -d 1 -x 1 -z 10 -j 5 -F $prefix.frgStore -f -o $prefix.fgbStore -L ofglist -I $prefix.ovlStore")) {
        rename "$prefix.cgb", "$prefix.cgb.FAILED";
        die "Failed.";
    }
    unlink "ofglist";
}


if (! -e "$prefix.cgi") {
    if (runCommand("$abin/consensus -P -U $prefix.frgStore $prefix.cgb")) {
        rename "$prefix.cgi", "$prefix.cgi.FAILED";
        die "Failed.\n";
    }
}

if (! -e "$prefix.cgw") {
    #  -CGW_CHECKPOINTS
    #  -j UNITIG_A_STATS (-j 1 from ian's script)
    if (runCommand("$abin/cgw -k 5 -r 4 -j 1 -s 2 -w 0 -T -P -f $prefix.frgStore -g $prefix.gkpStore -o $prefix $prefix.cgi")) {
        rename "$prefix.cgw", "$prefix.cgw.FAILED";
        die "Failed.\n";
    }
}

if (! -e "$prefix.cns") {
    system("cat $prefix.cgw $prefix.cgw_contigs $prefix.cgw_scaffolds > $prefix.cgw_total");
    if (runCommand("$abin/consensus -P $prefix.frgStore $prefix.cgw_total")) {
        rename "$prefix.cns", "$prefix.cns.FAILED";
        die "Failed.\n";
    }
    unlink "$prefix.cgw_total";
}


if (! -e "$prefix.asm") {
    if (runCommand("$abin/terminator -P -g $prefix.gkpStore -f $prefix.frgStore -i $prefix.cns -o $prefix.asm -m $prefix.map")) {
        rename "$prefix.asm", "$prefix.asm.FAILED";
        die "Failed\n";
    }
}

if (! -e "$prefix.fasta") {
    if (runCommand("$abin/process_scaffolds -f $prefix.scaffolds.fasta < $prefix.cns")) {
        rename "$prefix.fasta", "$prefix.fasta.FAILED";
        die "Failed.\n";
    }
    #if (runCommand("leaff -F $prefix.scaffolds.fasta -ii | sort -k4nr | head")) {
    #}
}






#  Utility to run a command and check the exit status
#
sub runCommand {
    my $cmd = shift @_;

    print "STARTING A COMMAND!-------------------------------------------------------------\n";
    system("date");
    print "$cmd\n";

    my $rc = 0xffff & system($cmd);

    system("date");
    print "FINISHED A COMMAND!-------------------------------------------------------------\n";

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!
    #
    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    }

    my $error = "ERROR: $cmd\n        failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } elsif ($rc > 0x80) {
        $rc >>= 8;
        $error .= "exit status $rc\n";
    } else {
        if ($rc & 0x80) {
            $rc &= ~0x80;
            $error .= "coredump from ";
        }
        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc]\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}
