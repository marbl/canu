#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

#
#  Mon Apr 22 11:29:23 EDT 2002 -- fixed asmversion/asmtaxon: C4 was hardcoded in
#  the xml generator.
#

$| = 1;

use strict;

#  We use the libBri package of usefull stuff.  It's located in the same place
#  as the ESTmapper.pl script.  That's what FindBin tells us.
#
use FindBin;
use lib "$FindBin::Bin";
use libBri;

my %p;
my @polishes;
my $lastdbID  = "NOT_A_VALID_DBID";
my $lastdbSCF = "NOT_A_VALID_SCAF";

my %ids;
my %map;

my $tier = "";
my $dir  = "";
my $uid  = 0;
my $tag  = "";

my $verbose = 0;

my $asmversion = "";
my $asmtaxon   = "";


############################################################
#
#  Useful subroutines
#


#  Returns a scaffold id given a defline.  If the defline is from a
#  Celera scaffold, it returns the proper thing.  Otherwise, it returns
#  a UID in the namespace provided on the command line
#
sub getSCFid {
    my $d = shift @_;
    if ($d =~ m/ga_uid=(\d+)\s/) {
        return("CELERA:$1");
    } elsif ($d =~ m/^>CRA\|(\d+)\s/) {
        return("CELERA:$1");
    } else {
        return(generateUID());
    }
}

#  Returns a cDNA id given a defline.  If the defline is from a
#  Celera sequence, it returns the proper thing.  Otherwise, it returns
#  a UID in the namespace provided on the command line
#
sub getESTid {
    my $d = shift @_;
    if ($d =~ m/^>CRA\|(\d+)\s/) {
        return("CELERA:$1");
    } else {
        return(generateUID());
    }
}

#  Returns a legal XML string, replacing the illegal characters (or,
#  some of them) with escaped versions.
#
sub XMLify {
    my $d = shift @_;
    $d =~ s/\&/&amp;/sg;
    $d =~ s/\"/&quot;/sg;
    $d =~ s/\'/&apos;/sg;
    $d =~ s/\>/&gt;/sg;
    $d =~ s/\</&lt;/sg;

    return($d);
}

#  Returns a UID in our namespace
#
sub generateUID {
    $uid++;
    return "$tag:$uid";
}

#  The usage
#
sub usage {
    print STDERR "usage: $0 [-verbose] -tier GBtier -tag UIDtag [-uid StartingUID] -dir outputDirectory\n";
    print STDERR "    -asmversion version -asmtaxon taxon\n";
    print STDERR "\n";
    print STDERR "The following aliases exist:\n";
    print STDERR "    -C3   == -asmversion  969031580 -asmtaxon \"Homo sapiens\"   (human     (C3) assembly)\n";
    print STDERR "    -C4   == -asmversion  985037830 -asmtaxon \"Homo sapiens\"   (human R26 (C4) assembly)\n";
    print STDERR "    -R27  == -asmversion 1006873711 -asmtaxon \"Homo sapiens\"   (human R27 (WGA5) assembly)\n";
    print STDERR "    -WGA2 == -asmversion  985359889 -asmtaxon \"Mus musculus\"   (mouse R12 (WGA2) assembly)\n";
    print STDERR "    -WGA3 == -asmversion  998319478 -asmtaxon \"Mus musculus\"   (mouse R13 (WGA3) assembly)\n";
    print STDERR "\n";
}


############################################################
#
#  Process the command line, checking for required arguments
#
while (scalar @ARGV > 0) {
    my $arg = shift @ARGV;

    if      ($arg eq "-tier") {
        $tier = shift @ARGV;
    } elsif ($arg eq "-tag") {
        $tag = shift @ARGV;
    } elsif ($arg eq "-uid") {
        $uid = shift @ARGV;
    } elsif ($arg eq "-dir") {
        $dir = shift @ARGV;
    } elsif ($arg eq "-verbose") {
        $verbose = 1;
    } elsif ($arg eq "-asmversion") {
        $asmversion = shift @ARGV;
    } elsif ($arg eq "-asmtaxon") {
        $asmtaxon = shift @ARGV;
    } elsif ($arg eq "-C3") {
        $asmversion = "969031580";
        $asmtaxon   = "Homo sapiens";
    } elsif ($arg eq "-C4") {
        $asmversion = "985037830";
        $asmtaxon   = "Homo sapiens";
    } elsif ($arg eq "-R27") {
        $asmversion = "1006873711";
        $asmtaxon   = "Homo sapiens";
    } elsif ($arg eq "-WGA2") {
        $asmversion = "985359889";
        $asmtaxon   = "Mus musculus";
    } elsif ($arg eq "-WGA3") {
        $asmversion = "998319478";
        $asmtaxon   = "Mus musculus";
    }
}

my $err = 0;

if ($tier eq "") {
    usage() if ($err == 0);
    print "ERROR: No tier given.\n";
    $err = 1;
}

if ($tag eq "") {
    usage() if ($err == 0);
    print "ERROR: No UID tag given.\n";
    $err = 1;
}

if ($dir eq "") {
    usage() if ($err == 0);
    print "ERROR: No output directory given.\n";
    $err = 1;
}

if ($asmversion eq "") {
    usage() if ($err == 0);
    print "ERROR: No assembly version given.\n";
    $err = 1;
}

if ($asmtaxon eq "") {
    usage() if ($err == 0);
    print "ERROR: No assembly taxon given.\n";
    $err = 1;
}

exit(1) if ($err);



############################################################
#
#  Most of the work is done in this function.  The driver is below it.
#
sub dumpXML {
    print "Found ", scalar @polishes, " polishes for scaffold $lastdbID / $lastdbSCF\n" if ($verbose);

    #  Hack the name
    #
    my $filename  = $lastdbSCF;
    my $subdir    = "";
    if ($filename =~ m/\D(\d+)(\d\d)$/) {

        #  Use the tier as a prefix on the filename, but we need
        #  to make it safe.  Non alphanumeric are removed.
        #
        $filename  =  "$tag:$1$2";
        $filename  =~ tr/[a-zA-Z0-9_+=\-:]//cds;
        $subdir    = "$2";
    } else {
        print STDERR "ERROR:  lastdbSCF format invalid -- it doesn't end in numbers!\n";
        print STDERR "ERROR:  lastdbSCF = '$filename'\n";
        return;
    }

    system("mkdir $dir/$subdir") if (! -d "$dir/$subdir");

    print "-> $dir/$subdir/$filename.gbf\n" if ($verbose);

    open(F, "> $dir/$subdir/$filename.gbf");

    print F "<game version=\"1.002\" assembly_version=\"$asmversion\" taxon=\"$asmtaxon\">\n";
    print F "<computational_analysis id=\"" . generateUID() . "\">\n";

    my ($day, $month, $year) = (localtime)[3,4,5];
    $month++;
    $year += 1900;

    print F " <date day=\"$day\" month=\"$month\" year=\"$year\"></date>\n";
    print F " <program>$tier</program>\n";
    print F " <type>Sim4_hit</type>\n";

    foreach my $P (@polishes) {
        print F " <result_set id=\"" . generateUID() . "\">\n";
        print F " <description>", XMLify($P->{'estDefLine'}), "</description>\n";

        my $strand = 1;
        if ($P->{'strandPrediction'} eq "unknown") {
            $strand = -1 if ($P->{'matchOrientation'} eq "complement");
        } else {
            $strand = -1 if ($P->{'strandPrediction'} eq "reverse");
        }

        foreach my $E (@{$P->{'exons'}}) {
            print F " <result_span id=\"" . generateUID() . "\">\n";
            print F "  <span_type>Sim4 Feature Detail</span_type>\n";

            my ($genB, $genE, $estB, $estE);

            #
            #  Convert all the coordinated to 0-space based
            #
            $genB = int($P->{'dbLo'} + $E->{'GENOMICstart'}) - 1;
            $genE = int($P->{'dbLo'} + $E->{'GENOMICend'});
            $estB = int($E->{'cDNAstart'}) - 1;
            $estE = int($E->{'cDNAend'});

            if ($strand == -1) {
                $genB = int($P->{'dbLo'} + $E->{'GENOMICend'});
                $genE = int($P->{'dbLo'} + $E->{'GENOMICstart'}) - 1;
                $estB = int($P->{'estLen'}) - int($E->{'cDNAend'});
                $estE = int($P->{'estLen'}) - int($E->{'cDNAstart'} - 1);
            }

            #if ($strand == -1) {
            #    $genB = int($P->{'dbLo'} + $E->{'GENOMICend'});
            #    $genE = int($P->{'dbLo'} + $E->{'GENOMICstart'}) - 1;
            #    $estB = int($P->{'estLen'}) + 1 - int($E->{'cDNAend'}) - 1;
            #    $estE = int($P->{'estLen'}) + 1 - int($E->{'cDNAstart'});
            #} else {
            #    $genB = int($P->{'dbLo'} + $E->{'GENOMICstart'}) - 1;
            #    $genE = int($P->{'dbLo'} + $E->{'GENOMICend'}) - 1;
            #    $estB = int($E->{'cDNAstart'});
            #    $estE = int($E->{'cDNAend'});
            #}

            #  The genomic description
            print F "  <seq_relationship id=\"$ids{$P->{'dbDefLine'}}\" type=\"query\">\n";
            print F "   <span><start>$genB</start><end>$genE</end></span>\n";

            if (defined($E->{'GENOMICalign'})) {
                print F "   <alignment>\n";
                if ($strand == 1) {
                    print F "   $E->{'GENOMICalign'}\n";
                } else {
                    my $x = reverse $E->{'GENOMICalign'};
                    $x =~ tr/acgtACGT/tgcaTGCA/;
                    print F "   $x\n";
                }
                print F "   </alignment>\n";
            }
            print F "  </seq_relationship>\n";

            #  The cDNA description
            print F "  <seq_relationship id=\"$ids{$P->{'estDefLine'}}\" type=\"sbjct\">\n";
            print F "   <span><start>$estB</start><end>$estE</end></span>\n";

            if (defined($E->{'cDNAalign'})) {
                print F "   <alignment>\n";
                if ($strand == 1) {
                    print F "   $E->{'cDNAalign'}\n";
                } else {
                    my $x = reverse $E->{'cDNAalign'};
                    $x =~ tr/acgtACGT/tgcaTGCA/;
                    print F "   $x\n";
                }
                print F "   </alignment>\n";
            }
            print F "  </seq_relationship>\n";

            my $alignLength = length($E->{'cDNAalign'});

            $_ = $E->{'cDNAalign'} . $E->{'GENOMICalign'};
            my $numGaps     = tr/-//;

            print F "  <property name=\"alignment_length\"   value=\"$alignLength\"></property>\n";
            print F "  <property name=\"num_identical\"      value=\"$E->{'numMatches'}\"></property>\n";
            print F "  <property name=\"num_similar\"        value=\"$E->{'numMatches'}\"></property>\n";
            print F "  <property name=\"num_gaps\"           value=\"$numGaps\"></property>\n";
            #print F "  <property name=\"score\"              value=\" ??? \"></property>\n";
            print F " </result_span>\n";
        }

        print F " <property name=\"subject_seq_id\"     value=\"$ids{$P->{'estDefLine'}}\"></property>\n";
        print F " <property name=\"subject_seq_length\" value=\"$P->{'estLen'}\"></property>\n";
        print F " </result_set>\n";
    }

    #print F "  <property name=\"score\"              value=\"\"></property>\n";
    print F "</computational_analysis>\n";
    print F "</game>\n";

    close(F);
}


############################################################
#
#  Begin of the main.
#
system("mkdir $dir") if (! -d $dir);

print STDERR "Converting.  Reading polishes sorted by scaffold from standard input.\n";
while (!eof(STDIN)) {
    undef %p;
    %p = &libBri::readPolish(*STDIN);

    if (defined($p{'raw'})) {
        if (int($p{'dbID'}) < ($lastdbID)) {
            print STDERR "$0: ERROR: The input polishes are not sorted by the scaffold id!\n";
            exit;
        }

        if ($p{'dbID'} ne $lastdbID) {
            if (scalar @polishes > 0) {
                dumpXML();
                print "Last UID used: '$tag:$uid'\n" if ($verbose);
            }
            
            undef %ids;
            undef %map;
            undef @polishes;

            $lastdbID  = $p{'dbID'};
            $lastdbSCF = getSCFid($p{'dbDefLine'});

            $ids{$p{'dbDefLine'}} = getSCFid($p{'dbDefLine'});
            $map{$p{'dbDefLine'}} = $p{'dbDefLine'};
        }

        if (!defined($ids{$p{'estDefLine'}})) {
            $ids{$p{'estDefLine'}} = getESTid($p{'estDefLine'});
            $map{$p{'estDefLine'}} = $p{'estDefLine'};
        }

        my $P = {};
        %{$P} = %p;
        push @polishes, $P;
    }
}

#  There is always something to dump, unless the input is empty,
#  in which case the user deserves bizarre output.
#
if (scalar @polishes > 0) {
    dumpXML();
}

print "Last UID used: '$tag:$uid'\n";
