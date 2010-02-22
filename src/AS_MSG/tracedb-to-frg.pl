#!/usr/bin/perl

#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2007, J. Craig Venter Instititue.
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
#  Takes as input a list of NCBI TraceDB formatted fasta, qual and xml
#  files.  Writes celera-assembler frg format (version 2 by default)
#  output.  Designed for converting LARGE amounts of data, using SGE.
#  Probably not very friendly.
#
#  Three steps:
#
#  1)  (parallel)   parse the xml.prefix.### files.
#  2)  (sequential) scan the parsed xml to find libraries and mates.
#                   Writes CA LIB and LKG messages.
#  3)  (parallel)   scan the fasta/qual files, write FRG messages.
#
#  If you have SGE installed, you can run all three stages, parallel if
#  appropriate, with:
#    tracedb-to-frg.pl -sgeoptions '-P project -A account etc' -sge xml*
#
#  Otherwise, you'll need to run (sequentially) with:
#    tracedb-to-frg.pl -xml xml  (once per xml file)
#    tracedb-to-frg.pl -lib xml* (ALL xml files)
#    tracedb-to-frg.pl -frg xml  (once per xml file)
#
#  Output goes in the CURRENT DIRECTORY.
#
#  Script diagnostic output from steps 1 and 3 should be empty; output
#  from step 2 will contain the mate rates.
#
#  The output is:
#
#  prefix.1.lib.frg -- all library info
#  prefix.2.###.frg -- fragment data, one per input fasta/qual pair
#  prefix.3.lkg.frg -- mate info
#
#
#  Other output:
#    *.lib -- a map from TA frag id to TA library id, along with some
#             info about the library
#    *.frglib -- a map from TA frag id to CA library UID
#
#

use strict;

#  Some centers decide to use very fine-grained SEQ_LIB_ID identifies
#  (the "ligationID"), which we use to group reads into a library.
#
#  This map will convert a SEQ_LIB_ID to a "library name".  It should be
#  many-to-one.  For example:
#    $seqLibIDFix{"829124933"} = "WUGSCsm1";
#    $seqLibIDFix{"829124934"} = "WUGSCsm1";
#
my %seqLibIDFix;

#  Switch deciding if we should use the SEQ_LIB_ID (1) or the
#  LIBRARY_ID (0) to get the name of the library.  Some centers choose
#  to not supply a SEQ_LIB_ID (the drosophila fragments from
#  Agencourt).
#
#  0 -- use LIBRARY_ID.
#  1 -- use SEQ_LIB_ID, but die when it isn't present.
#  2 -- use SEQ_LIB_ID, allow missing libs.
#
my $useSLI = 1;

#  This allows one to include only certain libraries in the output, or
#  to ignore specific libraries.  We'll do all the work, and just skip
#  the output steps.
#
my %seqLibExclude;
my %seqLibInclude;

#  Insetad of using the prefix of the ncbi files, we can define the
#  output name.  Do not define xmlNam, xmlIdx or tmpDir.
#
my $xmlNam;
my $xmlIdx;
my $tmpDir;
my $outNam;  # = "my-output-prefix";

my $version = 2;

my @sgefiles;
my $sgeOptions = undef;
my $xmlfile    = undef;
my @libfiles;
my $frgfile    = undef;
my $forNewbler = 0;
my $onlyOption = undef;

while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if      ($arg eq "-sge") {
        @sgefiles = @ARGV;
        undef @ARGV;
    } elsif ($arg eq "-sgeoptions") {
        $sgeOptions = shift @ARGV;
    } elsif ($arg eq "-xml") {
        $xmlfile = shift @ARGV;
    } elsif ($arg eq "-lib") {
        @libfiles = @ARGV;
        undef @ARGV;
    } elsif ($arg eq "-frg") {
        $forNewbler = 0;
        $frgfile = shift @ARGV;
    } elsif ($arg eq "-newbler") {
        $forNewbler = 1;
        $frgfile = shift @ARGV;
    } elsif ($arg eq "-only") {
        $outNam    = shift @ARGV;
        $seqLibInclude{$outNam} = 1;
        $onlyOption = "-only $outNam";
    } else {
        die "unknown option '$arg'.\n";
    }
}
if (!defined(@sgefiles) && !defined($xmlfile) && !defined(@libfiles) && !defined($frgfile)) {
    print STDERR "usage: $0 [-sge xml*]                          //  All steps on SGE\n";
    print STDERR "       $0 [-xml xml*] [-lib xml*] [-frg xml*]  //  Each stage independently\n";
    print STDERR "\n";
    print STDERR "       -sgeoptions \"options\"   supply this string to qsub\n";
    print STDERR "                                 (do NOT set -N; this will break the script)\n";
    print STDERR "       -newbler                  write format suitable for newbler\n";
    print STDERR "       -only LIBNAME             only process reads from library LIBNAME\n";
    print STDERR "                                 (may be used multiple times)\n";
    exit(1);
}

my $slexc = scalar(keys %seqLibExclude);
my $slinc = scalar(keys %seqLibInclude);

my $tcat = "";
#$tcat = "/home/bri/bin/tcat.osx |"             if (-x "/home/bri/bin/tcat.osx");
#$tcat = "/home/bri/bin/osFreeBSD-amd64/tcat |" if (-x "/home/bri/bin/osFreeBSD-amd64/tcat");
#$tcat = "/home/bri/bin/tcat.i386 |"            if (-x "/home/bri/bin/tcat.i386");

runSGE(@sgefiles) if (defined(@sgefiles));
runXML($xmlfile)  if (defined($xmlfile));
runLIB(@libfiles) if (defined(@libfiles));
runFRG($frgfile)  if (defined($frgfile) && ($forNewbler == 0));
runNBL($frgfile)  if (defined($frgfile) && ($forNewbler == 1));

#  everybody should have exited already.  if not, it's an error.
die "Someone didn't exit properly.\n";


#  These need to be global, they hold the ID of the next fasta/qual
#  sequence we're going to be reading, so that we can, ummm, well,
#  know what the next fasta/qual sequence that we're going to read is.
#
my $fhdr;
my $qhdr;

sub readXML () {
    my ($xid, $type, $template, $end, $lib, $libsize, $libstddev) = (".", ".", ".", ".", ".", ".", ".");
    my ($clrl, $clrr, $clvl, $clvr, $clql, $clqr, $clr, $clv, $clq);  #  need to be = undef!
    my $rgi;

    #  Early Baylor SeaUrchin data has more than two fragments for
    #  each SEQ_LIB_ID, so we also include the RUN_GROUP_ID, which
    #  seems to differentiate them.
    #
    #  Except it also breaks lots of mates in later files.

    $_ = <X>;
    while ($_ !~ m/^\s*<\/TRACE/i) {
        if (m!^\s*<TI>(\S+)<!i) {
            $xid = $1;
        }
        if (m!^\s*<TRACE_TYPE_CODE>(.*)<!i) {
            $type = $1;
            $type =~ tr/a-z/A-Z/;
        }
        if (m!^\s*<TEMPLATE_ID>(.*)</!i) {
            $template = $1;
        }
        if (m!^\s*<TRACE_END>(.*)<!i) {
            $end = $1;
        }
        if (m!^\s*<SEQ_LIB_ID>(.*)</!i) {
            if ($useSLI != 0) {
                $lib = $1;
                $lib = $seqLibIDFix{$lib} if (exists($seqLibIDFix{$lib}));
            }
        }
        if (m!^\s*<LIBRARY_ID>(.*)</!i) {
            if ($useSLI == 0) {
                $lib = $1;
            }
            if (($useSLI != 0) && ($lib eq ".")) {
                $lib = $1;
            }
        }

        #if (m!^\s*<RUN_GROUP_ID>(.*)</!i) {
        #    $rgi = $1;
        #}
        if (m!^\s*<INSERT_SIZE>(.*)</!i) {
            $libsize = $1;
        }
        if (m!^\s*<INSERT_STDEV>(.*)<!i) {
            $libstddev = $1;
        }

        if (m!^\s*<CLIP_LEFT>(.*)<!i) {
            $clrl = $1 - 1;
            $clrl = 0 if ($clrl < 0);
        }
        if (m!^\s*<CLIP_RIGHT>(.*)<!i) {
            $clrr = $1;
        }

        if (m!^\s*<CLIP_VECTOR_LEFT>(.*)<!i) {
            $clvl = $1 - 1;
            $clvl = 0 if ($clvl < 0);
        }
        if (m!^\s*<CLIP_VECTOR_RIGHT>(.*)<!i) {
            $clvr = $1;
        }

        if (m!^\s*<CLIP_QUALITY_LEFT>(.*)<!i) {
            $clql = $1 - 1;
            $clql = 0 if ($clql < 0);
        }
        if (m!^\s*<CLIP_QUALITY_RIGHT>(.*)<!i) {
            $clqr = $1;
        }

        $_= <X>;
    }

    #  We check for incomplete clear ranges later.
    $clr = "$clrl,$clrr";
    $clv = "$clvl,$clvr";
    $clq = "$clql,$clqr";

    $template =~ s/\s+/_/g;
    $lib      =~ s/\s+/_/g;

    $xid      =~ s/,//g;
    $lib      =~ s/,//g;

    return($xid, $type, $template, $end, $lib, $libsize, $libstddev, $clr, $clv, $clq);
}


sub readFasta ($) {
    my $convert = shift @_;
    my $fstr;

    while (!eof(F)) {
        $_ = <F>;
        chomp;

        if (m/^>/) {
            my $ret = $fhdr;

            if      (m/ti\|(\S+)\s/) {
                $fhdr = $1;
                $fhdr =~ s/,//g;
            } else {
                die "Failed to parse an ID out of the sequence defline '$_'\n";
            }

            if (defined($ret)) {
                return($ret, $fstr);
            }
        } else {
            my $q = $_;
            if ($convert) {
                $q =~ s/^\s+//;
                $q =~ s/\s+$//;
                $fstr .= $q;
            } else {
                $fstr .= "$q\n";
            }
        }
    }

    #   Got eof, return whatever
    my $ret = $fhdr;
    undef $fhdr;
    return($ret, $fstr);
}


sub readQual ($) {
    my $convert = shift @_;
    my $qstr;

    while (!eof(Q)) {
        $_ = <Q>;
        chomp;

        if (m/^>/) {
            my $ret = $qhdr;

            if      (m/ti\|(\S+)\s/) {
                $qhdr = $1;
                $qhdr =~ s/,//g;
            } else {
                die "Failed to parse an ID out of the quality defline '$_'\n";
            }

            if (defined($ret)) {
                return($ret, $qstr);
            }
        } else {
            my $q = $_;

            if ($convert) {
                $q =~ s/^\s+//;
                $q =~ s/\s+$//;
                foreach my $qv (split '\s+', $q) {
                    if ($qv > 60) {$qv = 60;}
                    $qstr .= chr(ord('0') + $qv);
                }
            } else {
                $qstr .= "$q\n";
            }
        }
    }

    #   Got eof, return whatever
    my $ret = $qhdr;
    undef $qhdr;
    return($ret, $qstr);
}




sub runSGE (@) {
    my @sgefiles = @_;

    foreach my $xmlfile (@sgefiles) {
        $xmlNam = $1 if ($xmlfile =~ m/xml.(.*)\.\d+.bz2$/);
        $xmlNam = $1 if ($xmlfile =~ m/xml.(.*)\.\d+.gz$/);
        $xmlNam = $1 if ($xmlfile =~ m/xml.(.*)\.\d+$/);
    }

    $tmpDir = "tafrg-$xmlNam";
    system("mkdir $tmpDir") if (! -e "$tmpDir");

    die "runSGE()-- Didn't find prefix.\n" if (!defined($xmlNam));

    my $jid = $xmlNam;
    if (defined($outNam)) {
        $jid = $outNam;
    }

    open(X, "> $tmpDir/tafrgX-$jid.sh") or die;
    open(L, "> $tmpDir/tafrgL-$jid.sh") or die;
    open(F, "> $tmpDir/tafrgF-$jid.sh") or die;

    print X "#!/bin/sh\n";
    print X "#\n";
    print X "#\$ -cwd\n";
    print X "#\$ -A assembly\n";
    print X "#\$ -j y\n";
    print X "#\n";

    print L "#!/bin/sh\n";
    print L "#\n";
    print L "#\$ -cwd\n";
    print L "#\$ -A assembly\n";
    print L "#\$ -j y\n";
    print L "#\n";
    print L "perl $0 $onlyOption -lib ";

    print F "#!/bin/sh\n";
    print F "#\n";
    print F "#\$ -cwd\n";
    print F "#\$ -A assembly\n";
    print F "#\$ -j y\n";
    print F "#\n";

    my $nf  = 0;

    foreach my $xmlfile (@sgefiles) {
        $nf++;

        print X "if [ \$SGE_TASK_ID -eq $nf ] ; then\n";
        print X "  perl $0 $onlyOption -xml $xmlfile\n";
        print X "fi\n";

        print L "\\\n  $xmlfile";

        print F "if [ \$SGE_TASK_ID -eq $nf ] ; then\n";
        print F "  perl $0 $onlyOption -frg $xmlfile\n";
        print F "fi\n";
    }

    print L "\n";

    close(X);
    close(L);
    close(F);

    system("chmod +x $tmpDir/tafrgX-$jid.sh");
    system("chmod +x $tmpDir/tafrgL-$jid.sh");
    system("chmod +x $tmpDir/tafrgF-$jid.sh");

    #  -b n  is the default; TIGR IT decided to change the default.

    system("qsub -b n -p -5 -N tfX-$jid                    -t 1-$nf -o $tmpDir/tafrgX-$jid-\\\$TASK_ID.out $sgeOptions ./$tmpDir/tafrgX-$jid.sh") and die;
    system("qsub -b n -p -4 -N tfL-$jid -hold_jid tfX-$jid          -o $tmpDir/tafrgL-$jid.out             $sgeOptions ./$tmpDir/tafrgL-$jid.sh") and die;
    system("qsub -b n -p -6 -N tfF-$jid -hold_jid tfL-$jid -t 1-$nf -o $tmpDir/tafrgF-$jid-\\\$TASK_ID.out $sgeOptions ./$tmpDir/tafrgF-$jid.sh") and die;

    exit(0);
}


sub runXML ($) {
    my $xmlfile = shift @_;
    my %types;

    if ($xmlfile =~ m/xml.(.*)\.(\d+).bz2$/) {
        $xmlNam = $1;
        $xmlIdx = $2;
        open(X, "bzip2 -dc xml.$xmlNam.$xmlIdx.bz2 | $tcat") or die "Failed to open xml.$xmlNam.$xmlIdx.bz2\n";
    }

    if ($xmlfile =~ m/xml.(.*)\.(\d+).gz$/) {
        $xmlNam = $1;
        $xmlIdx = $2;
        open(X, "gzip  -dc xml.$xmlNam.$xmlIdx.gz  | $tcat") or die "Failed to open xml.$xmlNam.$xmlIdx.gz\n";
    }

    if ($xmlfile =~ m/xml.(.*)\.(\d+)$/) {
        $xmlNam = $1;
        $xmlIdx = $2;
        open(X, "< xml.$xmlNam.$xmlIdx") or die "Failed to open xml.$xmlNam.$xmlIdx\n";
    }

    $tmpDir = "tafrg-$xmlNam";
    system("mkdir $tmpDir") if (! -e "$tmpDir");

    die "Didn't understand file '$xmlfile'.\n" if (!defined($xmlNam));

    if (-e "$tmpDir/$xmlNam.$xmlIdx.lib") {
        exit(0);
    }

    open(L, "> $tmpDir/$xmlNam.$xmlIdx.lib.tmp") or die "Failed to open '$tmpDir/$xmlNam.$xmlIdx.lib.tmp' for write.\n";

    while (!eof(X)) {
        my ($xid, $type, $template, $end, $lib, $libsize, $libstddev, $clr, $clv, $clq) = readXML();

        $types{$type}++;

        if (($type eq "WGS") ||
            ($type eq "SHOTGUN") ||
            ($type eq "CLONEEND") ||
            ($type eq "454")) {

            if ($lib eq ".") {
                print STDERR "No LIB found for $xid.\n";
                print STDERR "Consider setting useSLI.\n";
            }

            print L "$xid\t$template\t$end\t$lib\t$libsize\t$libstddev\t$clr\t$clv\t$clq\n";
        }
    }

    close(L);
    close(X);

    rename "$tmpDir/$xmlNam.$xmlIdx.lib.tmp", "$tmpDir/$xmlNam.$xmlIdx.lib";

    #foreach my $k (sort keys %types) {
    #    print "$xmlNam.$xmlIdx\t$k\t$types{$k}\n";
    #}

    exit(0);
}


sub runLIB (@) {
    my @libfiles = @_;

    my %mean;
    my %stddev;
    my %link;
    my %dist;
    my %distNeeded;

    foreach my $xmlfile (@libfiles) {

        if ($xmlfile =~ m/xml.(.*)\.(\d+).bz2$/) {
            $xmlNam = $1;
            $xmlIdx = $2;
        }

        if ($xmlfile =~ m/xml.(.*)\.(\d+).gz$/) {
            $xmlNam = $1;
            $xmlIdx = $2;
        }

        if ($xmlfile =~ m/xml.(.*)\.(\d+)$/) {
            $xmlNam = $1;
            $xmlIdx = $2;
        }

        die "runLIB()-- Didn't find the prefix in any of the input files (did you give me the xml?).\n" if (!defined($xmlNam));

        $tmpDir = "tafrg-$xmlNam";
        system("mkdir $tmpDir") if (! -e "$tmpDir");

        $outNam = $xmlNam if (!defined($outNam));

        if ((-e "$outNam.1.lib.frg") && (-e "$outNam.3.lkg.frg")) {
            print STDERR "Nothing to do.  '$outNam.1.lib.frg' and '$outNam.3.lkg.frg' exist.\n";
            print STDERR "(If later steps fail, remove these two files to recompute intermediate results)\n";
            exit(0);
        }

        open(F, "< $tmpDir/$xmlNam.$xmlIdx.lib") or die "Failed to open '$tmpDir/$xmlNam.$xmlIdx.lib'\n";

        while (<F>) {
            chomp;
            my ($id, $template, $end, $lib, $mn, $sd, $clr, $clv, $clq) = split '\t', $_;

            #  Keep track of all the libraries we've seen.
            $dist{$lib} = $lib;

            #  Skip it?
            next if (($slexc > 0) &&  exists($seqLibExclude{$lib}));
            next if (($slinc > 0) && !exists($seqLibInclude{$lib}));

            #  Keep track of the libraries we need to output
            $distNeeded{$lib}++;

            if (($mn < 1000) || ($sd == 0)) {
                undef $mn;
                undef $sd;
            }

            if (defined($mean{$lib}) && defined($mn)) {
                if ($mean{$lib} < $mn) {
                    $mean{$lib} = $mn;
                    $stddev{$lib} = $sd;
                }
            }

            if (!defined($mean{$lib}) && defined($mn)) {
                $mean{$lib} = $mn;
                $stddev{$lib} = $sd;
            }
        }

        close(F);
    }

    open(LIB, "> $outNam.1.lib.frg.tmp") or die "Failed to open LIB output '$outNam.1.lib.frg.tmp'\n";
    open(LKG, "> $outNam.3.lkg.frg.tmp") or die "Failed to open LKG output '$outNam.3.lkg.frg.tmp'\n";

    if ($version > 1) {
        print LIB "{VER\n";
        print LIB "ver:$version\n";
        print LIB "}\n";

        print LKG "{VER\n";
        print LKG "ver:$version\n";
        print LKG "}\n";
    }

    foreach my $distUID (keys %distNeeded) {

        #  Make a new bogus distance if the input never defined a valid one.
        #
        if (!defined($mean{$distUID})) {
            print STDERR "WARNING: creating bogus distance for library '$distUID'\n";
            $mean{$distUID}   = 2000;
            $stddev{$distUID} = 666;
        }

        #  And output.
        #
        if ($version == 1) {
            print LIB "{DST\n";
            print LIB "act:A\n";
            print LIB "acc:$distUID\n";
            print LIB "mea:$mean{$distUID}\n";
            print LIB "std:$stddev{$distUID}\n";
            print LIB "}\n";
        } else {
            print LIB "{LIB\n";
            print LIB "act:A\n";
            print LIB "acc:$distUID\n";
            print LIB "ori:I\n";
            print LIB "mea:$mean{$distUID}\n";
            print LIB "std:$stddev{$distUID}\n";
            print LIB "src:\n";
            print LIB ".\n";
            print LIB "nft:0\n";
            print LIB "fea:\n";
            print LIB ".\n";
            print LIB "}\n";
        }
    }

    foreach my $xmlfile (@libfiles) {
        my $numMate = 0;
        my $numFrag = 0;

        if ($xmlfile =~ m/xml.(.*)\.(\d+).bz2$/) {
            $xmlNam = $1;
            $xmlIdx = $2;
        }

        if ($xmlfile =~ m/xml.(.*)\.(\d+).gz$/) {
            $xmlNam = $1;
            $xmlIdx = $2;
        }

        if ($xmlfile =~ m/xml.(.*)\.(\d+)$/) {
            $xmlNam = $1;
            $xmlIdx = $2;
        }

        my $writeFrgLib = 0;
        if (! -e "$tmpDir/$xmlNam.$xmlIdx.frglib") {
            $writeFrgLib = 1;
            open(L, "> $tmpDir/$xmlNam.$xmlIdx.frglib.tmp") or die "Failed to open '$tmpDir/$xmlNam.$xmlIdx.frglib.tmp' for write.\n";
        }

        open(F, "< $tmpDir/$xmlNam.$xmlIdx.lib") or die "Failed to open '$tmpDir/$xmlNam.$xmlIdx.lib' for read\n";

        while (<F>) {
            my ($id, $template, $end, $lib, $mn, $sd, $clr, $clv, $clq) = split '\t', $_;

            #  Remember which library this frag came from.
            #
            if ($writeFrgLib) {
                if (! exists($dist{$lib})) {
                    print STDERR "WARNING: $lib doesn't exist in %dist.\n";
                }
                print L "$id\t$dist{$lib}\n";
            }

            #  Skip it?
            next if (($slexc > 0) &&  exists($seqLibExclude{$lib}));
            next if (($slinc > 0) && !exists($seqLibInclude{$lib}));

            #  Baylor sea urchin data, early, has more than one mate per
            #  template.  We assume those are well ordered, and just let
            #  them slip by as new mates.  The template originally (but
            #  has since been disabled in write-frag.pl) included the
            #  template and the run_group_id, until the run_group_id
            #  started changing across mates.
            #
            #  NOTE!  You can't use the extra sanity checking disabled
            #  below!
            #
            #($template, undef) = split '-', $template;

            if ($template ne ".") {
                if (defined($link{$template})) {

                    #  F -> 5'
                    #  R -> 3'

                    my $id1 = $id;
                    my ($id2, $end2) = split '\0', $link{$template};

                    if ($end eq "R") {
                        ($id1, $id2) = ($id2, $id1);
                    }

                    #  Check disabled because of Baylor sea urchin.
                    #if ($end2 eq $end) {
                    #    print STDERR "WARNING: template=$template has two $end reads.\n";
                    #}

                    if ($version == 1) {
                        print LKG "{LKG\n";
                        print LKG "act:A\n";
                        print LKG "typ:M\n";
                        print LKG "fg1:$id1\n";
                        print LKG "fg2:$id2\n";
                        print LKG "etm:0\n";
                        print LKG "dst:$dist{$lib}\n";
                        print LKG "ori:I\n";
                        print LKG "}\n";
                    } else {
                        print LKG "{LKG\n";
                        print LKG "act:A\n";
                        print LKG "frg:$id1\n";
                        print LKG "frg:$id2\n";
                        print LKG "}\n";
                    }

                    $numMate++;

                    delete $link{$template}
                } else {
                    $link{$template} = "$id\0$end";
                }
            }

            $numFrag++;
        }

        if ($writeFrgLib) {
            close(L);
            rename "$tmpDir/$xmlNam.$xmlIdx.frglib.tmp", "$tmpDir/$xmlNam.$xmlIdx.frglib";
        }

        close(F);

        my $percMated = 0;
        if ($numFrag > 0) {
            $percMated = int(10000 * 2 * $numMate / $numFrag) / 100;
        }

        print STDERR "$xmlNam.$xmlIdx: frags=$numFrag links=$numMate ($percMated%)\n";
    }

    close(LKG);
    close(LIB);

    rename "$outNam.1.lib.frg.tmp", "$outNam.1.lib.frg";
    rename "$outNam.3.lkg.frg.tmp", "$outNam.3.lkg.frg";

    exit(0);
}


sub runFRG ($) {
    my $frgfile = shift @_;

    if ($frgfile =~ m/xml.(.*)\.(\d+).bz2$/) {
        $xmlNam = $1;
        $xmlIdx  = $2;
        open(X, "bzip2 -dc xml.$xmlNam.$xmlIdx.bz2   | $tcat") or die "Failed to open 'xml.$xmlNam.$xmlIdx.bz2' for read\n";
        open(F, "bzip2 -dc fasta.$xmlNam.$xmlIdx.bz2 | $tcat") or die "Failed to open 'fasta.$xmlNam.$xmlIdx.bz2' for read\n";
        open(Q, "bzip2 -dc qual.$xmlNam.$xmlIdx.bz2  | $tcat") or die "Failed to open 'qual.$xmlNam.$xmlIdx.bz2' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+).gz$/) {
        $xmlNam = $1;
        $xmlIdx  = $2;
        open(X, "gzip -dc xml.$xmlNam.$xmlIdx.gz   | $tcat") or die "Failed to open 'xml.$xmlNam.$xmlIdx.gz' for read\n";
        open(F, "gzip -dc fasta.$xmlNam.$xmlIdx.gz | $tcat") or die "Failed to open 'fasta.$xmlNam.$xmlIdx.gz' for read\n";
        open(Q, "gzip -dc qual.$xmlNam.$xmlIdx.gz  | $tcat") or die "Failed to open 'qual.$xmlNam.$xmlIdx.gz' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+)$/) {
        $xmlNam = $1;
        $xmlIdx  = $2;
        open(X, "< xml.$xmlNam.$xmlIdx")    or die "Failed to open 'xml.$xmlNam.$xmlIdx' for read\n";
        open(F, "< fasta.$xmlNam.$xmlIdx")  or die "Failed to open 'fasta.$xmlNam.$xmlIdx' for read\n";
        open(Q, "< qual.$xmlNam.$xmlIdx")   or die "Failed to open 'qual.$xmlNam.$xmlIdx' for read\n";
    }

    $tmpDir = "tafrg-$xmlNam";
    system("mkdir $tmpDir") if (! -e "$tmpDir");

    #print STDERR "lib == $tmpDir/$xmlNam.$xmlIdx.frglib\n";
    open(L, "< $tmpDir/$xmlNam.$xmlIdx.frglib") or die "Failed to open '$tmpDir/$xmlNam.$xmlIdx.frglib' for read\n";

    $outNam = $xmlNam if (!defined($outNam));

    if (-e "$outNam.2.$xmlIdx.frg.bz2") {
        exit(0);
    }

    #  Count how many fragments are in this lib.  If zero, just stop now.
    #
    {
        my $fragsHere = 0;

        open(LT, "< $tmpDir/$xmlNam.$xmlIdx.frglib") or die "Failed to open '$tmpDir/$xmlNam.$xmlIdx.frglib' for read\n";
        while (<LT>) {
            my $lid = <LT>;
            my $lib;

            ($lid, $lib) = split '\s+', $lid;

            #  Skip it?.
            #
            next if (($slexc > 0) &&  exists($seqLibExclude{$lib}));
            next if (($slinc > 0) && !exists($seqLibInclude{$lib}));

            $fragsHere++;
        }
        close(LT);

        if ($fragsHere == 0) {
            print STDERR "No fragments we care about found in $tmpDir/$xmlNam.$xmlIdx.frglib.\n";
            exit(0);
        }

        print STDERR "Found $fragsHere in $tmpDir/$xmlNam.$xmlIdx.frglib.\n";
    }

    #open(FRG, "| $tcat bzip2 -9c > $outNam.2.$xmlIdx.frg.bz2.tmp") or die "Failed to open '$outNam.2.$xmlIdx.frg.bz2.tmp' for write\n";
    open(FRG, "> $outNam.2.$xmlIdx.frg.tmp") or die "Failed to open '$outNam.2.$xmlIdx.frg.tmp' for write\n";

    my $numFrags = 0;
    my $haveMore = 1;
    $haveMore  = !eof(X);
    $haveMore &= !eof(F);
    $haveMore &= !eof(Q);
    $haveMore &= !eof(L);

    while ($haveMore) {
        my ($xid, $type, $template, $end, $lib, $libsize, $libstddev, $clr, $clv, $clq) = readXML();
        my ($sid, $seq) = readFasta(1);
        my ($qid, $qlt) = readQual(1);

        $haveMore  = !eof(X);
        $haveMore &= !eof(F);
        $haveMore &= !eof(Q);
        $haveMore &= !eof(L);

        next if (($type ne "WGS") &&
                 ($type ne "SHOTGUN") &&
                 ($type ne "CLONEEND") &&
                 ($type ne "454"));

        my $lid = <L>;
        ($lid, $lib) = split '\s+', $lid;

        #  Skip it?.
        #
        next if (($slexc > 0) &&  exists($seqLibExclude{$lib}));
        next if (($slinc > 0) && !exists($seqLibInclude{$lib}));

        my $len = length $seq;

        $clv = undef      if ($clv eq ",");
        $clq = undef      if ($clq eq ",");
        $clr = undef      if ($clr eq ",");

        #  Incomplete clear; WUGSC Dros Yakuba only gave the left
        #  vector clear.
        #
        $clv = "0$clv"    if ($clv =~ m/^,\d+$/);
        $clq = "0$clq"    if ($clq =~ m/^,\d+$/);
        $clr = "0$clr"    if ($clr =~ m/^,\d+$/);

        $clv = "$clv$len" if ($clv =~ m/^\d+,$/);
        $clq = "$clq$len" if ($clq =~ m/^\d+,$/);
        $clr = "$clr$len" if ($clr =~ m/^\d+,$/);

        #  If no clr defined, make up one from clq and clv.  Note that 
        #  clq is not propagated to the gkpStore.
        #
        if (!defined($clr)) {
            my $l = 0;
            my $r = $len;

            if (defined($clv)) {
                my ($a,$b) = split ',', $clv;
                $l = $a if ($l < $a);
                $r = $b if ($b < $r);
            }
            if (defined($clq)) {
                my ($a,$b) = split ',', $clq;
                $l = $a if ($l < $a);
                $r = $b if ($b < $r);
            }
            $clr = "$l,$r";
        }

        #  Because some centers apparently cannot count, we need
        #  to check the clear ranges.
        #
        if ($clv =~ m/^(\d+),(\d+)$/) { $clv = "$1,$len" if ($2 > $len); }
        if ($clq =~ m/^(\d+),(\d+)$/) { $clq = "$1,$len" if ($2 > $len); }
        if ($clr =~ m/^(\d+),(\d+)$/) { $clr = "$1,$len" if ($2 > $len); }

        if (($xid eq $sid) && ($xid eq $qid) && ($xid eq $lid)) {
            if ($version == 1) {
                print FRG "{FRG\n";
                print FRG "act:A\n";
                print FRG "acc:$xid\n";
                print FRG "typ:R\n";
                print FRG "src:\n.\n";
                print FRG "etm:0\n";
                print FRG "seq:\n$seq\n.\n";
                print FRG "qlt:\n$qlt\n.\n";
                print FRG "clr:$clr\n";
                print FRG "}\n";
            } else {
                if ($numFrags == 0) {
                    print FRG "{VER\n";
                    print FRG "ver:$version\n";
                    print FRG "}\n";
                }
                print FRG "{FRG\n";
                print FRG "act:A\n";
                print FRG "acc:$xid\n";
                print FRG "rnd:1\n";
                print FRG "sta:G\n";
                print FRG "lib:$lib\n";
                print FRG "pla:0\n";
                print FRG "loc:0\n";
                print FRG "src:\n.\n";
                print FRG "seq:\n$seq\n.\n";
                print FRG "qlt:\n$qlt\n.\n";
                print FRG "hps:\n.\n";
                print FRG "clv:$clv\n" if (defined($clv));
                #print FRG "clq:$clq\n" if (defined($clq));
                print FRG "clr:$clr\n";
                print FRG "}\n";
            }
            $numFrags++;
        } else {
            print STDERR "ID mismatch: X='$xid' =?= S='$sid'\n";
            print STDERR "ID mismatch: X='$xid' =?= Q='$qid'\n";
            print STDERR "ID mismatch: X='$xid' =?= L='$lid'\n";
            die;
        }
    }

    close(L);
    close(Q);
    close(F);
    close(X);
    close(FRG);

    #rename "$outNam.2.$xmlIdx.frg.bz2.tmp", "$outNam.2.$xmlIdx.frg.bz2";
    #unlink "$outNam.2.$xmlIdx.frg.bz2" if ($numFrags == 0);
    rename "$outNam.2.$xmlIdx.frg.tmp", "$outNam.2.$xmlIdx.frg";
    unlink "$outNam.2.$xmlIdx.frg" if ($numFrags == 0);

    exit(0);
}



sub runNBL ($) {
    my $frgfile = shift @_;

    if ($frgfile =~ m/xml.(.*)\.(\d+).bz2$/) {
        $xmlNam = $1;
        $xmlIdx  = $2;
        open(X, "bzip2 -dc xml.$xmlNam.$xmlIdx.bz2   | $tcat") or die "Failed to open 'xml.$xmlNam.$xmlIdx.bz2' for read\n";
        open(F, "bzip2 -dc fasta.$xmlNam.$xmlIdx.bz2 | $tcat") or die "Failed to open 'fasta.$xmlNam.$xmlIdx.bz2' for read\n";
        open(Q, "bzip2 -dc qual.$xmlNam.$xmlIdx.bz2  | $tcat") or die "Failed to open 'qual.$xmlNam.$xmlIdx.bz2' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+).gz$/) {
        $xmlNam = $1;
        $xmlIdx  = $2;
        open(X, "gzip -dc xml.$xmlNam.$xmlIdx.gz   | $tcat") or die "Failed to open 'xml.$xmlNam.$xmlIdx.gz' for read\n";
        open(F, "gzip -dc fasta.$xmlNam.$xmlIdx.gz | $tcat") or die "Failed to open 'fasta.$xmlNam.$xmlIdx.gz' for read\n";
        open(Q, "gzip -dc qual.$xmlNam.$xmlIdx.gz  | $tcat") or die "Failed to open 'qual.$xmlNam.$xmlIdx.gz' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+)$/) {
        $xmlNam = $1;
        $xmlIdx  = $2;
        open(X, "< xml.$xmlNam.$xmlIdx")    or die "Failed to open 'xml.$xmlNam.$xmlIdx' for read\n";
        open(F, "< fasta.$xmlNam.$xmlIdx")  or die "Failed to open 'fasta.$xmlNam.$xmlIdx' for read\n";
        open(Q, "< qual.$xmlNam.$xmlIdx")   or die "Failed to open 'qual.$xmlNam.$xmlIdx' for read\n";
    }

    $tmpDir = "tafrg-$xmlNam";
    system("mkdir $tmpDir") if (! -e "$tmpDir");

    open(L, "< $tmpDir/$xmlNam.$xmlIdx.frglib") or die "Failed to open '$tmpDir/$xmlNam.$xmlIdx.frglib' for read\n";

    $outNam = $xmlNam if (!defined($outNam));

    if (-e "$outNam.2.$xmlIdx.frg.bz2") {
        exit(0);
    }

    open(FNA, "> $outNam.2.$xmlIdx.fna.tmp")      or die "Failed to open '$outNam.2.$xmlIdx.fna.tmp' for write\n";
    open(QNA, "> $outNam.2.$xmlIdx.fna.qual.tmp") or die "Failed to open '$outNam.2.$xmlIdx.fna.qual.tmp' for write\n";

    my $haveMore = 1;
    $haveMore  = !eof(X);
    $haveMore &= !eof(F);
    $haveMore &= !eof(Q);
    $haveMore &= !eof(L);

    while ($haveMore) {
        my ($xid, $type, $template, $end, $lib, $libsize, $libstddev, $clr, $clv, $clq) = readXML();
        my ($sid, $seq) = readFasta(0);
        my ($qid, $qlt) = readQual(0);

        $end = "F" if ($end eq "FORWARD");
        $end = "R" if ($end eq "REVERSE");

        if (($type eq "WGS") ||
            ($type eq "SHOTGUN") ||
            ($type eq "CLONEEND") ||
            ($type eq "454")) {
            my $lll = <L>;
            my $lid;
            ($lid, $lib) = split '\s+', $lll;

            my $len = length $seq;

            #  Incomplete clear; WUGSC Dros Yakuba only gave the left
            #  vector clear.
            #
            $clv = "0$clv"    if ($clv =~ m/^,\d+$/);
            $clq = "0$clq"    if ($clq =~ m/^,\d+$/);
            $clr = "0$clr"    if ($clr =~ m/^,\d+$/);

            $clv = "$clv$len" if ($clv =~ m/^\d+,$/);
            $clq = "$clq$len" if ($clq =~ m/^\d+,$/);
            $clr = "$clr$len" if ($clr =~ m/^\d+,$/);

            $clr = "0,$len" if (!defined($clr));

            #  Because some centers apparently cannot count, we need
            #  to check the clear ranges.
            #
            if ($clv =~ m/^(\d+),(\d+)$/) { $clv = "$1,$len" if ($2 > $len); }
            if ($clq =~ m/^(\d+),(\d+)$/) { $clq = "$1,$len" if ($2 > $len); }
            if ($clr =~ m/^(\d+),(\d+)$/) { $clr = "$1,$len" if ($2 > $len); }

            $clv =~ s/,/-/;
            $clq =~ s/,/-/;
            $clr =~ s/,/-/;

            #  Newbler is picky about clear ranges.  It wants them to
            #  actually include sequence.  This ugly construction
            #  omits reads with empty clear ranges.
            #
            if ($clr =~ m/(\d+)-(\d+)/) {

                #  Newbler uses base-based clear ranges; this script
                #  uses space-based clear.
                #
                my $cll = $1 + 1;
                my $clr = $2;

                if ($1 < $2) {
                    if (($xid eq $sid) && ($xid eq $qid) && ($xid eq $lid)) {
                        print FNA ">$xid template=$template dir=$end library=$lib trim=$clr\n$seq";
                        print QNA ">$xid template=$template dir=$end library=$lib trim=$clr\n$qlt";
                    } else {
                        print STDERR "ID mismatch: X='$xid' =?= S='$sid'\n";
                        print STDERR "ID mismatch: X='$xid' =?= Q='$qid'\n";
                        print STDERR "ID mismatch: X='$xid' =?= L='$lid'\n";
                        die;
                    }
                }
            }
        }

        $haveMore  = !eof(X);
        $haveMore &= !eof(F);
        $haveMore &= !eof(Q);
        $haveMore &= !eof(L);
    }

    close(L);
    close(Q);
    close(F);
    close(X);
    close(FNA);
    close(QNA);

    rename "$outNam.2.$xmlIdx.fna.tmp",      "$outNam.2.$xmlIdx.fna";
    rename "$outNam.2.$xmlIdx.fna.qual.tmp", "$outNam.2.$xmlIdx.fna.qual";

    exit(0);
}
