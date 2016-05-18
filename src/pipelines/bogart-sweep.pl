#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2016-MAR-10
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

my $wrk = "/work/canuassemblies/sent";
my $asm = "test";

system("mkdir -p $wrk")  if (! -d $wrk);

my $gs = "5000000";
my $b  = 6000;

my $m  = 4;
my $t  = 1;

my $d = "all";

my (@EG, @EB, @EM, @ER, @OL, @RS, @NS, @CS);

@EG = (  "0.0100",  "0.0200",  "0.0300",  "0.0400",  "0.0500",  "0.0600",  "0.0700",  "0.0800",  "0.0900",  "0.1000" );
@EB = ( "-0.0050", "+0.0050" );
@EM = ( "-0.0050", "+0.0050" );
@ER = ( "-0.0050", "+0.0050" );
@OL = ( "100", "500", "1000", "2500", "5000", "7500", "10000" );
@RS = ( "-no", "-RS" );
@NS = ( "-no", "-NS" );
@CS = ( "-no", "-CS" );

@EG = (  "0.0400",  "0.0500",  "0.0600" );
@EB = ( "+0.0125" );
@EM = ( "-0.0125" );
@ER = ( "+0.0000" );
@OL = ( "50", "500", "1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000", "11000" );
@RS = ( "-RS" );
@NS = ( "-NS" );
@CS = ( "-CS" );

@EG = (  "0.0500" );
@EB = ( "+0.0125" );
@EM = ( "-0.0125" );
@ER = ( "+0.0000" );
@OL = ( "4100", "4200", "4300", "4400", "4500", "4600", "4700", "4800", "4900", "5100", "5200", "5300", "5400", "5500", "5600", "5700", "5800", "5900" );
@RS = ( "-RS" );
@NS = ( "-NS" );
@CS = ( "-CS" );

undef @OL;
for (my $ii=1100; $ii<1200; $ii += 1) {
    push @OL, $ii;
}

@OL = ( "1135", "1136", "1137", "1138" );


#-unassembled 2 1000 0.75 0.75 2 -repeatdetect 6 11 15 -threads 1 -D most

foreach my $eg (@EG) {
foreach my $eb (@EB) {
foreach my $em (@EM) {
foreach my $er (@ER) {
foreach my $ol (@OL) {
foreach my $rs (@RS) {
foreach my $ns (@NS) {
foreach my $cs (@CS) {
    my ($egl, $ebl, $eml, $erl, $oll) = ($eg, $eb, $em, $er, $ol);

    $ebl = $egl + $1  if ($eb =~ m/^\+(\d+.\d+)/);
    $ebl = $egl - $1  if ($eb =~ m/^-(\d+.\d+)/);

    $eml = $egl + $1  if ($em =~ m/^\+(\d+.\d+)/);
    $eml = $egl - $1  if ($em =~ m/^-(\d+.\d+)/);

    $erl = $egl + $1  if ($er =~ m/^\+(\d+.\d+)/);
    $erl = $egl - $1  if ($er =~ m/^-(\d+.\d+)/);

    $egl = sprintf("%6.4f", $egl);
    $ebl = sprintf("%6.4f", $ebl);
    $eml = sprintf("%6.4f", $eml);
    $erl = sprintf("%6.4f", $erl);

    $oll = sprintf("%05d", $ol);

    my $path = "test-eg$egl-eb$ebl-em$eml-er$erl-ol$oll$rs$ns$cs";

    print "$path\n";

    system("mkdir -p $path")  if (! -d $path);

    open(F, "> $wrk/$path/bogart.sh") or die "can't open '$wrk/$path/bogart.sh' for writing: $!\n";
    print F "#!/bin/sh\n";
    print F "\n";
    print F "cd $wrk/$path\n";
    print F "\n";
    print F "if [ ! -e test.tigStore ] ; then\n";
    print F "  /work/canu/FreeBSD-amd64/bin/bogart \\\n";
    print F "    -G $wrk/$asm.gkpStore \\\n";
    print F "    -O $wrk/$asm.ovlStore \\\n";
    print F "    -T test.tigStore -o test\\\n";
    print F "    -B $b -M $m -threads $t \\\n";
    print F "    -gs $gs \\\n";
    print F "    -eg $egl -eb $ebl -em $eml -er $erl -el $ol \\\n";
    print F "    -RS \\\n"  if ($rs eq "-RS");
    print F "    -NS \\\n"  if ($ns eq "-NS");
    print F "    -CS \\\n"  if ($cs eq "-CS");
    print F "    -unassembled 2 1000 0.75 0.75 2 \\\n";
    print F "    -repeatdetect 6 32 5 \\\n";
    print F "    -D $d \\\n"  if (length($d) > 0);
    print F "  > bogart.err 2>& 1\n";
    print F "fi\n";
    print F "\n";
    close(F);

    open(F, "> $wrk/$path/utgcns.sh") or die "can't open '$wrk/$path/utgcns.sh' for writing: $!\n";
    print F "#!/bin/sh\n";
    print F "\n";
    print F "cd $wrk/$path\n";
    print F "\n";
    print F "if [ ! -e test.fasta ] ; then\n";
    print F "  /work/canu/FreeBSD-amd64/bin/utgcns \\\n";
    print F "    -G $wrk/$asm.gkpStore \\\n";
    print F "    -T test.tigStore 1 . \\\n";
    print F "    -O test.cns -L test.lay -A test.fasta\n";
    print F "fi\n";
    print F "\n";
    print F "rm -f test.tigStore/seqDB.v002.dat\n";
    print F "rm -f test.tigStore/seqDB.v002.tig\n";
    print F "\n";
    print F "/work/canu/FreeBSD-amd64/bin/tgStoreLoad \\\n";
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  -T test.tigStore 2 \\\n";
    print F "  test.cns\n";
    print F "\n";
    print F "/work/canu/FreeBSD-amd64/bin/tgStoreDump \\\n";
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  -T test.tigStore 2 \\\n";
    print F "  -consensus -fasta -contigs -bubbles \\\n";
    print F "> contigs.fasta\n";
    print F "\n";
    print F "rm -f *.delta\n";
    print F "rm -f *.coords\n";
    print F "rm -f *.png\n";
    print F "\n";
    print F "sh /work/scripts/dotplot.sh usmarc /data/references/salmonella_enterica_usmarc_3124.1-cp006631.1.fasta contigs.fasta\n";
    print F "sh /work/scripts/dotplot.sh serge  /data/references/salmonella_enterica_serge.fasta contigs.fasta\n";
    print F "\n";
    print F "cp -fp usmarc.png $wrk/$path.usmarc.png\n";
    print F "cp -fp serge.png  $wrk/$path.serge.png\n";
    print F "\n";
    close(F);

    system("sh $wrk/$path/bogart.sh");
    system("qsub -q vomit.q -cwd -j y -o /dev/null $wrk/$path/utgcns.sh > /dev/null 2>&1");
}
}
}
}
}
}
}
}

