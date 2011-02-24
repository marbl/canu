#!/usr/bin/perl

use strict;

my $gkpStore;
my $tigStore;

my $ctgStoreVersion = 2;
my $utgStoreVersion = 2;

my %deletedFragments;


sub deleteUnitig($) {
    my $utgID    = shift @_;

    my $saveDir  = "deleted-unitig-$utgID";

    system("mkdir $saveDir") if (! -e $saveDir);

    system("tigStore -g $gkpStore -t $tigStore $utgStoreVersion -u $utgID -d layout          > $saveDir/$utgID.layout");
    system("tigStore -g $gkpStore -t $tigStore $utgStoreVersion -u $utgID -d consensus       > $saveDir/$utgID.consensus");
    system("tigStore -g $gkpStore -t $tigStore $utgStoreVersion -u $utgID -d consensusgappes > $saveDir/$utgID.consensusgappes");
    system("tigStore -g $gkpStore -t $tigStore $utgStoreVersion -u $utgID -d multialign      > $saveDir/$utgID.multialign");

    open(F, "> $saveDir/$utgID.layout.delete");
    print F "unitig $utgID\n";
    print F "len 0\n";
    print F "cns\n";
    print F "qlt\n";
    print F "data.unitig_coverage_stat 0.0\n";
    print F "data.unitig_microhet_prob 1.0\n";
    print F "data.unitig_status        X\n";
    print F "data.unitig_unique_rept   X\n";
    print F "data.contig_status        U\n";
    print F "data.num_frags            0\n";
    print F "data.num_unitigs          0\n";
    close(F);

    system("tigStore -g $gkpStore -t $tigStore $utgStoreVersion -R $saveDir/$utgID.layout.delete");

    open(F, "< $saveDir/$utgID.layout");
    while (<F>) {
        if (m/^FRG\s+type\s+.\s+ident\s+(\d+)\s+/) {
            $deletedFragments{$1}++;
        }
    }
    close(F);
}



sub deleteContig($) {
    my $ctgID    = shift @_;

    my $saveDir  = "deleted-contig-$ctgID";

    system("mkdir $saveDir") if (! -e $saveDir);

    system("tigStore -g $gkpStore -t $tigStore $ctgStoreVersion -c $ctgID -d layout          > $saveDir/$ctgID.layout");
    system("tigStore -g $gkpStore -t $tigStore $ctgStoreVersion -c $ctgID -d consensus       > $saveDir/$ctgID.consensus");
    system("tigStore -g $gkpStore -t $tigStore $ctgStoreVersion -c $ctgID -d consensusgappes > $saveDir/$ctgID.consensusgappes");
    system("tigStore -g $gkpStore -t $tigStore $ctgStoreVersion -c $ctgID -d multialign      > $saveDir/$ctgID.multialign");

    open(F, "> $saveDir/$ctgID.layout.delete");
    print F "contig $ctgID\n";
    print F "len 0\n";
    print F "cns\n";
    print F "qlt\n";
    print F "data.unitig_coverage_stat 0.0\n";
    print F "data.unitig_microhet_prob 1.0\n";
    print F "data.unitig_status        X\n";
    print F "data.unitig_unique_rept   X\n";
    print F "data.contig_status        U\n";
    print F "data.num_frags            0\n";
    print F "data.num_unitigs          0\n";
    close(F);

    system("tigStore -g $gkpStore -t $tigStore $ctgStoreVersion -R $saveDir/$ctgID.layout.delete");

    open(F, "< $saveDir/$ctgID.layout");
    while (<F>) {
        if (m/unitig\s+(\d+)/) {
            deleteUnitig($1);
        }
        if (m/^FRG\s+type\s+.\s+ident\s+(\d+)\s+/) {
            $deletedFragments{$1}++;
        }
    }
    close(F);
}



sub deleteFragments() {

    #  Remove the mate association from the fragments to delete.  This has the side
    #  effect of printing the mate IID.

    open(F, "> deleted-fragments.1.getmate");
    foreach my $frg (sort { $a <=> $b } keys %deletedFragments) {
        print F "frg iid $frg mateiid 0\n";
    }
    close(F);

    system("gatekeeper --edit deleted-fragments.1.getmate $gkpStore > deleted-fragments.1.getmate.out");

    #  Read the mate IIDs, and remove the mate from those too.

    open(I, "< deleted-fragments.1.getmate.out");
    open(F, "> deleted-fragments.2.delmate");
    while (<I>) {
        if (m/frg\s+uid\s+.*\s+mateiid\s+(\d+)\s+->\s+mateiid\s+\d+/) {
            print F "frg iid $1 mateiid 0\n";
        }
    }
    close(I);
    close(F);

    system("gatekeeper --edit deleted-fragments.2.delmate $gkpStore > deleted-fragments.2.delmate.out");

    #  Then, finally, delete the fragments.

    open(F, "> deleted-fragments.3.delfrag");
    foreach my $frg (sort { $a <=> $b } keys %deletedFragments) {
        print F "frg iid $frg isdeleted 1\n";
    }
    close(F);

    system("gatekeeper --edit deleted-fragments.3.delfrag $gkpStore > deleted-fragments.3.delfrag.out");
}



my $utgID;
my $ctgID;

my $err      = 0;

while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if      ($arg eq "-g") {
        $gkpStore = shift @ARGV;

    } elsif ($arg eq "-t") {
        $tigStore = shift @ARGV;

    } elsif ($arg eq "-u") {
        $utgID = shift @ARGV;

    } elsif ($arg eq "-c") {
        $ctgID = shift @ARGV;

    } else {
        $err++;
        print STDERR "Invalid arg '$arg'\n";
    }
}

if (($err > 0) ||
    (!defined($gkpStore)) ||
    (!defined($tigStore)) ||
    ($ctgID + $utgID == 0)) {
    print STDERR "usage: $0 -g gkpStore -t tigStore [-c ctgID | -u utgID]\n";
    print STDERR "\n";
    print STDERR "  Deletes a contig (and all associated unitigs) or a single unitig from\n";
    print STDERR "  the given tigStore.  This operation will invalidate ALL SCAFFOLDING.\n";
    print STDERR "\n";
    print STDERR "  The deleted objects and other temporary files are saved in directories\n";
    print STDERR "  named:\n";
    print STDERR "    deleted-contig-(ID)\n";
    print STDERR "    deleted-unitig-(ID)\n";
    print STDERR "  Deleting a contig also delete the unitigs it contains, with each unitig\n";
    print STDERR "  saved in a separate directory.\n";
    exit(1);
}

deleteContig($ctgID)  if ($ctgID > 0);
deleteUnitig($utgID)  if ($utgID > 0);

deleteFragments();
