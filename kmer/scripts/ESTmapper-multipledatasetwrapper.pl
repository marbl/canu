#!/usr/local/bin/perl

$| = 1;

use strict;
use FindBin;
use lib "$FindBin::Bin";
use libBri;

my $ests = "/dev5/walenz/svi-human-input/map_ests_orphans.fasta";
my $mrna = "/dev5/walenz/svi-human-input/map_mrna_orphans.fasta";


#
#  Configure the list of things to map
#


sub configureThingsToMap {
    my $p = {};

    $p->{'ests'}   = shift;
    $p->{'mrna'}   = shift;
    $p->{'genome'} = shift;
    $p->{'label'}  = shift;
    $p->{'mask'}   = shift;
    
    return(@_, $p);
}


sub runMapper {
    my $p = shift;
    my $masking;

    print "------------------------------------------------------------\n";

    print "ESTs: '$p->{'label'}'\n" if (defined($p->{'ests'}));
    print "mRNA: '$p->{'label'}'\n" if (defined($p->{'mrna'}));

    $masking = "-maskmers \"\" ";
    $masking = "-maskmers $p->{'mask'}.fasta " if (defined($p->{'mask'}));

    if (defined($p->{'mask'})) {
        if (! -e "$p->{'mask'}.fasta") {
            print "Building frequentMers in '$p->{'mask'}'\n";

            my $meryl = "/home/walenzbp/projects/meryl/meryl";

            system("$meryl -v -B -f -m 20 -t 28 -s $p->{'genome'} -o $p->{'mask'}");
            system("$meryl -v -Dt -n 1000 -s $p->{'mask'} > $p->{'mask'}.fasta");
        } else {
            print "frequentMers found in '$p->{'mask'}'\n";
        }
    }

    if (defined($p->{'ests'})) {
        system("nice -20 /usr/local/bin/perl /prod/IR01/ESTmapper/ESTmapper.pl " .
               "-mapest " .
               "/dev5/walenz/svi-orphanests-mapped-to-$p->{'label'} " .
               "$p->{'ests'} " .
               "$p->{'genome'} " .
               "-verbose -stats -localsearches 4 -searchthreads 4 -localpolishes 4 $masking");
    }

    if (defined($p->{'mrna'})) {
        system("nice -20 /usr/local/bin/perl /prod/IR01/ESTmapper/ESTmapper.pl " .
               "-mapmrna " .
               "/dev5/walenz/svi-orphanmrna-mapped-to-$p->{'label'} " .
               "$p->{'mrna'} " .
               "$p->{'genome'} " .
               "-v -stats -localsearches 4 -searchthreads 4 -localpolishes 4 $masking");
    }
}





my @thingsToMap;



@thingsToMap = configureThingsToMap($ests,
                                    undef,
                                    "/dev5/walenz/GENOMES/CHGD_assembly_ge10k_R27.fasta",
                                    "R27",
                                    "/dev5/walenz/frequentMers-R27-20",
                                    @thingsToMap);



foreach my $p (@thingsToMap) {
    runMapper($p);
}


