#!/usr/bin/perl -w
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
#
# Usage:
#
# arg1 :
#
#######################################################################
#
# The modules of the Celera Assembler have a set of common command line
# arguments that are used in the production mode.  In fact, all
# command line arguments that are intended to be used in the production
# shown have unique meanings.
#
#######################################################################

use strict;
use Getopt::Long;
use Cwd;

# $prod = cwd();  # path to working directory

my $AS_ROOT = $ENV{AS_ROOT};
my $AS_BIN = $ENV{AS_BIN};
my $TIMER="/usr/bin/time";

# $AS_BIN = $ARGV[1];
# $prefix = $ARGV[0];

#my $AS_BIN = "";
my $prefix = "";
my $restart = 0;
my $run_consensus = 1;
my $graph_history = 1; # This records the major unitigger graph operations in the file name.

my $REPEAT_LIB="${AS_ROOT}/lib/$prefix.lib";
# if [ ! -e ${REPEAT_LIB} ]; then
#  my $REPEAT_LIB="${AS_ROOT}/lib/drosophila_repeats.lib";
# fi

my $OUTPUT_MODE = "-P";
my $AnalOpt = "-A 2";
my $survey_of_graphs = 0;

sub make_ovl_store
{
    my $OLD_OUTPUT_MODE = $OUTPUT_MODE;
    my $OLD_CREATE_MODE="-f";
    my $NEW_CREATE_MODE="-c -f";
    # The -f option is used to clobber an existant store with the same name.
    # The -c option is used to create an empty store.
    
    # my $CELSIM_CMD="${AS_BIN}/celsim";
    my $GATEKEEPER_CMD="${AS_BIN}/gatekeeper -X -C -N -Q ${OUTPUT_MODE} ${OLD_CREATE_MODE}";
    # my $URCSCREENER_CMD="${AS_BIN}/urc_screener -ra ${OLD_OUTPUT_MODE} ${REPEAT_LIB}";
    my $URCSCREENER_CMD="${AS_BIN}/urc_screener -r -f -s ${OLD_OUTPUT_MODE} ${REPEAT_LIB}";
    my $OVERLAP_CMD="${AS_BIN}/overlap -w ${OUTPUT_MODE} ${OLD_CREATE_MODE}";
    
    my $cmdline;
    
    if($restart <= 11) { 
        system("date");
        print "Remove bin.\n";
        $cmdline="rm -rf bin";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    
    if($restart <= 12) { 
        system("date");
        print "Make bin.\n";
        $cmdline="mkdir bin";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    if($restart <= 13) { 
        system("date");
        print "Fill bin.\n";
        $cmdline="cp $AS_ROOT/bin/* bin";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    if($restart <= 14) { 
        system("date");
        print "Fill bin.\n";
        $cmdline = "cp $AS_ROOT/scripts/make_partitionFile bin";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
        $cmdline = "cp $AS_ROOT/scripts/submit_consensus bin";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
        $AS_BIN = "./bin";
    }
    
    print "Began the Celera Assembler at `date` in `pwd`.\n";
    print "By the way, this script depends on AS_ROOT=$AS_ROOT\n\n";
    
    #$number_of_fragments="`grep "{FRG" $prefix.frg | wc -l `"
    #print "The input file $prefix.frg has $number_of_fragments fragment records.\n";
    
    if($restart <= 50) {
        print "Count the number of overlaps in the batch\n";
        $cmdline="countmessages < $prefix.frg > $prefix.cnt";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    
    if($restart <= 51) {
        my $number_of_fragments=`gawk '/FRG/ { print \$3}' $prefix.cnt`;
        print "The input file $prefix.frg has $number_of_fragments fragment records.\n";
    }
    
    if($restart <= 60) {
        print "\n=========== GKP: gatekeeper =====================\n";
        print "  $prefix.gkpStore is the gatekeeper store directory\n";
        $cmdline="$TIMER $GATEKEEPER_CMD $prefix.gkpStore $prefix.frg 2> $prefix.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    
    if($restart <= 70) {
        print "\n=========== URT: screen for ubiquitous repeats ==\n";
        if( 0 ) {
            print "  $REPEAT_LIB is the screen item library\n";
            $cmdline="$TIMER $URCSCREENER_CMD $prefix.inp 2>> $prefix.log";
        } else {
            print "  Skipping ....\n";
            $cmdline="ln -f -s $prefix.inp $prefix.urc";
        }
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    } 
    
    if($restart <= 80) {
        print "\n=========== OVL: Create overlap records =========\n";
        print "  $prefix.frgStore is the fragment store directory\n";
        $cmdline="$TIMER $OVERLAP_CMD $prefix.frgStore $prefix.urc 2>> $prefix.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    
    if($restart <= 90) {
        print "\n=========== OVL STORE: Create overlap store =========\n";
        $cmdline="$TIMER ${AS_BIN}/grow-olap-store -cf -o $prefix.ovlStore $prefix.ovl 2>> $prefix.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    
}


sub make_fgb_graphs
{
#    my $ReaperFileOptions = "-x 1 -z 10 -e 0.015";
    my $ReaperFileOptions = "-x 1 -z 10";
    my $ReaperOptions = "$ReaperFileOptions -d 0 -M 0";
    my $FGBOptions = "$ReaperFileOptions -d 1 -M 1";
    
    # Using append mode it would be:
    #$AS_BIN/make_OFG_from_FragStore $prefix.frgStore > $prefix.ofg;
    #$AS_BIN/fgb -P -A 1 -x 1 -z 10 -d 0 -M 0 -I dros.ovlStore -Q 0 -f -c -o $prefix.reaperStore $prefix.ofg;
    #$AS_BIN/fgb -P -A 1 -x 1 -z 10 -d 0 -M 0 -I dros.ovlStore -Q 1 -a -i $prefix.reaperStore;
    #$AS_BIN/cgb -P -A 1 -b 0 -j 5 $prefix.frgStore $prefix.reaperStore;


    my $cmdline;
    
    if($restart <= 110) { 
        system("date");
        print "Make OFG messages.\n";
        $cmdline = "$AS_BIN/make_OFG_from_FragStore $prefix.frgStore > $prefix.ofg";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 or die("Failed: $cmdline");
    }
    
    if(($restart <= 130)) {
        system("date");
        $cmdline = "$AS_BIN/fgb -P $AnalOpt $ReaperOptions -I $prefix.ovlStore -Q 0 -f -c -o ${prefix}_RQ0.fgbStore $prefix.ofg 1> ${prefix}_RQ0.fgbStore.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
    }
    
    if(($restart <= 140) && $survey_of_graphs) {
        system("date");
        $cmdline = "$AS_BIN/fgb -P $AnalOpt $ReaperOptions -I $prefix.ovlStore -Q 1 -i ${prefix}_RQ0.fgbStore -f -o ${prefix}_RQ0_RQ1.fgbStore 1> ${prefix}_RQ0_RQ1.fgbStore.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
    }
    
    if(($restart <= 150) && $survey_of_graphs) {
        system("date");
        $cmdline = "$AS_BIN/fgb -P $AnalOpt $ReaperOptions -I $prefix.ovlStore -Q 2 -i ${prefix}_RQ0.fgbStore -f -o ${prefix}_RQ0_RQ2.fgbStore 1> ${prefix}_RQ0_RQ2.fgbStore.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
    }
    
    if(($restart <= 160) && $survey_of_graphs) {
        system("date");
        $cmdline = "$AS_BIN/fgb -P $AnalOpt $ReaperOptions -I $prefix.ovlStore -Q 2 -i ${prefix}_RQ0_RQ1.fgbStore -f -o ${prefix}_RQ0_RQ1_RQ2.fgbStore 1> ${prefix}_RQ0_RQ1_RQ2.fgbStore.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
    }
    

    if(($restart <= 170)) {
        system("date");
        $cmdline = "$AS_BIN/fgb -P $AnalOpt $FGBOptions -Q 0 -i ${prefix}_RQ0.fgbStore -f -o ${prefix}_RQ0_T.fgbStore 1> ${prefix}_RQ0_T.fgbStore.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
    }

    if(($restart <= 180) && $survey_of_graphs) {
        system("date");
        $cmdline = "$AS_BIN/fgb -P $AnalOpt $FGBOptions -Q 0 -i ${prefix}_RQ0_RQ1.fgbStore -f -o ${prefix}_RQ0_RQ1_T.fgbStore 1> ${prefix}_RQ0_RQ1_T.fgbStore.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
    }

    print "The alternate graphs (before smoothing) are made.\n";
}


sub assemble_the_graphs
{
    # Global parameters: $restart, $prefix
    my $a_graph = $_[0];

    # Temporary variables:
    my $cmdline;
    my $outgraph = "${a_graph}"; # default result
    my $instore;
    
    my $CNSopts = "";
    #my $CGW_CMD="${AS_BIN}/cgw -c -j 1 -k 5 -w 1 -s 2 -r 4 -T ${OUTPUT_MODE}";
    my $CGW_CMD="${AS_BIN}/cgw -c -j 1 -k 5 -r 4 -s 2 -w 1 -T ${OUTPUT_MODE}";
    my $CONSENSUS_CMD="${AS_BIN}/consensus ${OUTPUT_MODE} ${CNSopts}";
    my $TERMINATOR_CMD="${AS_BIN}/terminator ${OUTPUT_MODE}";
    my $TESTER_CMD="${AS_BIN}/tester -h";

    my $smooth_spurs = 1;
    my $smooth_bubbles = 1;
    
    if( $smooth_spurs && $restart <= 210) {
        # Spur fragment smoothing
        system("date");
        $cmdline = "$AS_BIN/cgb -P $AnalOpt -s -b 0 -j 5 -U 0 $prefix.frgStore $a_graph.fgbStore 2>> $a_graph.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
        
        my $num_crappies = 0;
        
        system("date");
        $num_crappies = `gawk '{ if(\$0 !~ /^#/) print \$0;}' ${a_graph}.cgb_crappies | wc -l'`;
        print "$num_crappies crappies\n";
        
        if( $num_crappies > 0 ) {
            $cmdline = "$AS_BIN/repair_breakers -S ${a_graph}.cgb_crappies -F ${prefix}.frgStore -O ${a_graph}.breaker_overlaps.ovl ${a_graph}.cgb  2>> $a_graph.log";
            print "Executing: $cmdline\n";
            system($cmdline) == 0 || die("$cmdline");
            ###############################################
            
            ###############################################
            # fgb (2nd time)
            $instore = "${a_graph}.fgbStore";
            if( $graph_history ) { $outgraph = "${a_graph}_spur_T";}
            my $outstore = "${outgraph}.fgbStore";
            $cmdline="${AS_BIN}/fgb -P $AnalOpt -d 1 -M 1 -i $instore -f -o $outstore ${a_graph}.breaker_overlaps.ovl  2>> $a_graph.log";
            print "Executing: $cmdline\n";
            system($cmdline) == 0 || die("$cmdline");
            $a_graph = "${outgraph}";
            ###############################################
        }
    }

    if($smooth_bubbles && $restart <= 250) {
        # Bubble smoothing
        ###############################################
        # cgb (2nd time)
        $instore = "${a_graph}.fgbStore";
        $cmdline = "${AS_BIN}/cgb -P $AnalOpt -j 5 -b 0 -U 1 ${prefix}.frgStore $instore 2>> $a_graph.log";
        # run_cgb_with_bubbles.pl -C "$AnalOpt -j 5 -b 0" ${prefix}.frgStore ${prefix}.fgbStore
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
        
        ###############################################
        # fgb (2nd time)
        if( $graph_history ) { $outgraph = "${a_graph}_bubble_T";}
        my $outstore = "${outgraph}.fgbStore";
        $cmdline = "${AS_BIN}/fgb -P $AnalOpt -d 1 -M 1 -i $instore -f -o $outstore ${a_graph}.bubble_edges.ovl 2>> $a_graph.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
        $a_graph = "${outgraph}";
        ###############################################
    }
    
    
    if($restart <= 310) {
        system("date");
        $cmdline = "$AS_BIN/cgb -P $AnalOpt -b 0 -j 5 -U 0 $prefix.frgStore $a_graph.fgbStore  2>> $a_graph.log";
        print "Executing: $cmdline\n";
        system($cmdline) == 0 || die("$cmdline");
    }

    if($run_consensus){    
        if(1) {
#        if($restart <= 320) {
#            $cmdline = "ln -s $a_graph.cgb $prefix.cgb";
#            print "Executing: $cmdline\n";
#            system($cmdline) == 0 or die("Failed: $cmdline");
#        }
            if($restart <= 330) {
                system("date");
                print "Run consensus\n";
                $cmdline = "$AS_BIN/consensus -P -U $prefix.frgStore ${a_graph}.cgb  2>> $a_graph.log";
                print "Executing: $cmdline\n";
                system($cmdline) == 0 or die("Failed: $cmdline");
            }
            
        } else {
            if($restart <= 340) { 
                system("date");
                print "Make partition file\n";
                $cmdline = "$AS_BIN/make_partitionFile $prefix  2>> $a_graph.log";
                print "Executing: $cmdline\n";
                system($cmdline) == 0 || die("make_partitionFile died");
            }
            
            if($restart <= 350) { 
                system("date");
                print "Make partitioned fragment store\n";
                $cmdline = "$AS_BIN/partitionFragStore $prefix.partFile $prefix.frgStore $prefix.frgStore_part  2>> $a_graph.log";
                print "Executing: $cmdline\n";
                system($cmdline) == 0 || die("partitionFragStore");
            }
            
            if($restart <= 360) { 
                system("date");
                print "Submit consensus\n";
                $cmdline = "$AS_BIN/submit_consensus assembly $a_graph  2>> $a_graph.log";
                print "Executing: $cmdline\n";
                system($cmdline) == 0 or die("Failed: $cmdline");
            }
        }
        
        if($restart <= 410) {
            print "\n=========== CGW: scaffold and contig  =============\n";
            system("date");
            print "Run cgw\n";
            $cmdline = "$CGW_CMD -f $prefix.frgStore -g $prefix.gkpStore -o $a_graph $a_graph.cgi 2>> $a_graph.log";
            print "Executing: $cmdline\n";
            system($cmdline) == 0 or die("Failed: $cmdline");
        }
        
        if($restart <= 511) {
            print "\n=========== CNS: consensus builder  =============\n";
            system("date");
            $cmdline="$TIMER $CONSENSUS_CMD $prefix.frgStore $a_graph.cgw 2>> $a_graph.log";
            print "Executing: $cmdline\n";
            system($cmdline) == 0 or die("Failed: $cmdline");
        }
        
# print "\n=========== TST: tester =========================\n";
# $cmdline="$TIMER $TESTER_CMD $prefix 2>> $a_graph.log";
# print "Executing: $cmdline\n";
# system($cmdline) == 0 or die("Failed: $cmdline");
        
        if($restart <= 711) {
            print "\n=========== TER: terminator =====================\n";
            system("date");
            $cmdline="$TIMER $TERMINATOR_CMD -g $prefix.gkpStore -f $prefix.frgStore -i $a_graph.cns -o $a_graph.asm -m $a_graph.map 2>> $a_graph.log";
            print "Executing: $cmdline\n";
            system($cmdline) == 0 or die("Failed: $cmdline");
        }
        
        if(0) {
            if($restart <= 811) {
                system("date");
                print "Clean the directory of working files *.inp *.ovl *.cgb *cgw.\n";
                print "Assembler finished at `date`\n";
            }
        }
    } # if $run_consensus
    system("date");
}



sub usage
{
    print "Required arguments:\n";
    print "   -bin=<string>\n";
    print "   -prefix=<string>\n";
    print "Optional arguments:\n";
    print "   -restart=<integer>\n";
}

sub main
{
    
    &GetOptions(
                "bin=s", => \$AS_BIN, 
                "prefix=s", => \$prefix,
                "restart=i", => \$restart );

    print "AS_BIN   = $AS_BIN\n";
    print "prefix   = $prefix\n";
    print "restart  = $restart\n";

    if(
       ($AS_BIN eq "") ||
       ($prefix eq "") )
    {
        usage();
        exit 1;
    }

    make_ovl_store();
    make_fgb_graphs();

#    my @the_graph_list;
#    if( $survey_of_graphs ) {
#        @the_graph_list = ( "${prefix}_RQ0", "${prefix}_RQ0_RQ1", "${prefix}_RQ0_RQ2",
#                        "${prefix}_RQ0_T", "${prefix}_RQ0_RQ1_T" );
#    } elsif {
       my @the_graph_list = ( "${prefix}_RQ0_T" );
#    }
    
    my $a_graph;
    foreach $a_graph (@the_graph_list) {
        assemble_the_graphs($a_graph);
    }
}

main();
exit 0;

