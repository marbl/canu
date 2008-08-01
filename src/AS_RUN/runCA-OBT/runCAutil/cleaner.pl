use strict;

#  Assembly all done, remove some of the crud.

sub cleaner () {
    my $cleanType = getGlobal("cleanup");
    my $cleanValu = 0;

    print STDERR "The Cleaner has arrived.  Doing '$cleanType'.\n";

    $cleanValu = 0  if ($cleanType =~ m/none/);
    $cleanValu = 1  if ($cleanType =~ m/light/);
    $cleanValu = 2  if ($cleanType =~ m/heavy/);
    $cleanValu = 3  if ($cleanType =~ m/aggressive/);


    if ($cleanValu >= 1) {
        #
        #  Remove some of the more useless output files,
        #  and many of the stores and whatnot that can be recreated.
        #
        rmrf("$asm.obtStore");
        rmrf("0-mercounts/*blocks", "0-mercounts/*sequence");
        rmrf("0-overlaptrim-overlap/overlap*out");
        rmrf("1-overlapper/overlap*out");
        rmrf("4-unitigger/$asm.fge", "4-unitigger/$asm.fgv");
        rmrf("7*/rezlog");
    }


    if ($cleanValu >= 2) {
        #
        #
        #
    }


    if ($cleanValu >= 3) {
        #
        #  Nuke everything except 9-terminator.  Be paranoid about doing it.
        #
        rmrf("0-mercounts");
        rmrf("0-overlaptrim");
        rmrf("0-overlaptrim-overlap");
        rmrf("1-overlapper");
        rmrf("2-frgcorr");
        rmrf("3-ovlcorr");
        rmrf("4-unitigger");
        rmrf("5-consensus");
        rmrf("7-[0-9]-CGW");
        rmrf("7-[0-9]-ECR");
        rmrf("7-CGW");
        rmrf("8-consensus");
        rmrf("$asm.SeqStore");
        rmrf("$asm.asm");
        rmrf("$asm.frg");
        rmrf("$asm.gkpStore");
        rmrf("$asm.obtStore");
        rmrf("$asm.ovlStore");
        rmrf("$asm.qc");
    }


    if ($cleanType =~ m/compress/) {
        #  Compress *.err (usually tiny)
        #  Compress overlaps (*ovb)
        #  Compress checkpoints (*ckp.*[0-9])
    }
}
