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
#  appropriate with:
#    tracedb-to-frg.pl -sgeoptions '-P project -A account etc' -sge xml*
#
#  Otherwise, you'll need to run (sequentially) with:
#    tracedb-to-frg.pl -xml xml  (once per file)
#    tracedb-to-frg.pl -lib xml* (ALL xml files)
#    tracedb-to-frg.pl -frg xml  (once per file)
#
#  Output goes in the CURRENT DIRECTORY.
#
#  Script diagnostic output from steps 1 and 3 should be empty; output
#  from step 2 will contain the mate rates.
#
#  The output is:
#
#  hypocrea_jecorina.1.lib.frg -- all library info
#  hypocrea_jecorina.2.001.frg -- fragment data
#  hypocrea_jecorina.2.002.frg --   one per input fasta/qual pair
#  hypocrea_jecorina.3.lkg.frg -- mate info
#
#
#  Other output:
#    *.lib -- a map from TA frag id to TA library id, along with some
#             info about the library
#    *.frglib -- a map from TA frag id to CA library UID
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

my $version = 2;

my @sgefiles;
my $sgeOptions = undef;
my $xmlfile    = undef;
my @libfiles;
my $frgfile    = undef;
my $forNewbler = 0;

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
    } else {
        die "unknown option '$arg'.\n";
    }
}
if (!defined(@sgefiles) && !defined($xmlfile) && !defined(@libfiles) && !defined($frgfile)) {
    print STDERR "usage: $0 [-sge xml*]                          //  All steps on SGE\n";
    print STDERR "       $0 [-xml xml*] [-lib xml*] [-frg xml*]  //  Each stage independently\n";
    exit(1);
}

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
        }
        if (m!^\s*<TEMPLATE_ID>(.*)</!i) {
            $template = $1;
        }
        if (m!^\s*<TRACE_END>(.*)<!i) {
            $end = $1;
        }
        if (m!^\s*<SEQ_LIB_ID>(.*)</!i) {
            if ($lib eq ".") {
                $lib = $1;

                if (defined($seqLibIDFix{$lib})) {
                    $lib = $seqLibIDFix{$lib};
                }
            }
        }
        if (m!^\s*<LIBRARY_ID>(.*)</!i) {
            $lib = $1;
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

    #  Yes, OR.  We check for incomplete clear ranges later.
    $clr = "$clrl,$clrr" if (defined($clrl) || defined($clrr));
    $clv = "$clvl,$clvr" if (defined($clvl) || defined($clvr));
    $clq = "$clql,$clqr" if (defined($clql) || defined($clqr));

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

    my $prefix;
    foreach my $xmlfile (@sgefiles) {
        $prefix = $1 if ($xmlfile =~ m/xml.(.*)\.\d+.bz2$/);
        $prefix = $1 if ($xmlfile =~ m/xml.(.*)\.\d+.gz$/);
        $prefix = $1 if ($xmlfile =~ m/xml.(.*)\.\d+$/);
    }

    die "runSGE()-- Didn't find prefix.\n" if (!defined($prefix));

    open(X, "> tafrgX-$prefix.sh") or die;
    open(L, "> tafrgL-$prefix.sh") or die;
    open(F, "> tafrgF-$prefix.sh") or die;

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
    print L "perl $0 -lib ";

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
        print X "  perl $0 -xml $xmlfile\n";
        print X "fi\n";

        print L "\\\n  $xmlfile";

        print F "if [ \$SGE_TASK_ID -eq $nf ] ; then\n";
        print F "  perl $0 -frg $xmlfile\n";
        print F "fi\n";
    }

    print L "\n";

    close(X);
    close(L);
    close(F);

    system("chmod +x tafrgX-$prefix.sh");
    system("chmod +x tafrgL-$prefix.sh");
    system("chmod +x tafrgF-$prefix.sh");

    #  -b n  is the default; TIGR IT decided to change the default.

    system("qsub -b n -p -5 -N tafrgX-$prefix                          -t 1-$nf -o tafrgX-$prefix-\\\$TASK_ID.out $sgeOptions ./tafrgX-$prefix.sh") and die;
    system("qsub -b n -p -4 -N tafrgL-$prefix -hold_jid tafrgX-$prefix          -o tafrgL-$prefix.out             $sgeOptions ./tafrgL-$prefix.sh") and die;
    system("qsub -b n -p -6 -N tafrgF-$prefix -hold_jid tafrgL-$prefix -t 1-$nf -o tafrgF-$prefix-\\\$TASK_ID.out $sgeOptions ./tafrgF-$prefix.sh") and die;

    exit(0);
}


sub runXML ($) {
    my $xmlfile = shift @_;
    my $prefix;
    my %types;

    if ($xmlfile =~ m/xml.(.*\.\d+).bz2$/) {
        $prefix = $1;
        open(X, "bzip2 -dc xml.$prefix.bz2 | $tcat") or die "Failed to open xml.$prefix.bz2\n";
    }

    if ($xmlfile =~ m/xml.(.*\.\d+).gz$/) {
        $prefix = $1;
        open(X, "gzip  -dc xml.$prefix.gz  | $tcat") or die "Failed to open xml.$prefix.gz\n";
    }

    if ($xmlfile =~ m/xml.(.*\.\d+)$/) {
        $prefix = $1;
        open(X, "< xml.$prefix") or die "Failed to open xml.$prefix\n";
    }

    die "Didn't understand file '$xmlfile'.\n" if (!defined($prefix));

    if (-e "$prefix.lib") {
        exit(0);
    }

    open(L, "> $prefix.lib.tmp") or die "Failed to open '$prefix.lib.tmp' for write.\n";

    while (!eof(X)) {
        my ($xid, $type, $template, $end, $lib, $libsize, $libstddev, $clr, $clv, $clq) = readXML();

        $types{$type}++;

        if (($type eq "WGS") ||
            ($type eq "SHOTGUN") ||
            ($type eq "CLONEEND") ||
            ($type eq "454")) {
            print L "$xid\t$template\t$end\t$lib\t$libsize\t$libstddev\t$clr\t$clv\t$clq\n";
        }
    }

    close(L);
    close(X);

    rename "$prefix.lib.tmp", "$prefix.lib";

    foreach my $k (sort keys %types) {
        print "$prefix\t$k\t$types{$k}\n";
    }

    exit(0);
}


sub runLIB (@) {
    my @libfiles = @_;

    my %mean;
    my %stddev;
    my %link;
    my %dist;
    my %distNeeded;
    my $distUID = 6660000;  #  Starting UID for distances

    my $numDist = 0;
    my $numMate = 0;
    my $numFrag = 0;

    my $prefix;

    foreach my $xmlfile (@libfiles) {

        if ($xmlfile =~ m/xml.(.*\.\d+).bz2$/) {
            $prefix = $1;
        }

        if ($xmlfile =~ m/xml.(.*\.\d+).gz$/) {
            $prefix = $1;
        }

        if ($xmlfile =~ m/xml.(.*\.\d+)$/) {
            $prefix = $1;
        }

        die "runLIB()-- Didn't find the prefix in any of the input files (did you give me the xml?).\n" if (!defined($prefix));
        open(F, "< $prefix.lib")    or die "Failed to open '$prefix.lib'\n";

        while (<F>) {
            chomp;
            my ($id, $template, $end, $lib, $mn, $sd, $clr, $clv, $clq) = split '\s+', $_;

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
                $numDist++;
            }

            $distNeeded{$lib}++;
        }

        close(F);
    }

    if ($prefix =~ m/^(.*)\.\d+$/) {
        $prefix = $1;
    } else {
        die "Malformed prefix '$prefix'\n";
    }

    if ((-e "$prefix.1.lib.frg") && (-e "$prefix.3.lkg.frg")) {
        exit(0);
    }

    open(LIB, "> $prefix.1.lib.frg.tmp") or die "Failed to open LIB output '$prefix.1.lib.frg.tmp'\n";
    open(LKG, "> $prefix.3.lkg.frg.tmp") or die "Failed to open LKG output '$prefix.3.lkg.frg.tmp'\n";

    if ($version > 1) {
        print LIB "{VER\n";
        print LIB "ver:$version\n";
        print LIB "}\n";

        print LKG "{VER\n";
        print LKG "ver:$version\n";
        print LKG "}\n";
    }

    foreach my $d (keys %distNeeded) {

        #  Make a new bogus distance if the input never defined a valid one.
        #
        if (!defined($mean{$d})) {
            print STDERR "WARNING: creating bogus distance for library '$d'\n";
            $mean{$d} = 2000;
            $stddev{$d} = 666;
        }

        #  HACK!  Use the original library id for the distUID; otherwise,
        #  $distUID is a global and updated for each new distance record.
        #
        #$distUID = $d;
        
        #  And output.
        #
        if ($version == 1) {
            print LIB "{DST\n";
            print LIB "act:A\n";
            print LIB "acc:$distUID\n";
            print LIB "mea:$mean{$d}\n";
            print LIB "std:$stddev{$d}\n";
            print LIB "}\n";
        } else {
            print LIB "{LIB\n";
            print LIB "act:A\n";
            print LIB "acc:$distUID\n";
            print LIB "ori:I\n";
            print LIB "mea:$mean{$d}\n";
            print LIB "std:$stddev{$d}\n";
            print LIB "src:\n";
            print LIB "NCBI:$d\n";
            print LIB ".\n";
            print LIB "nft:0\n";
            print LIB "fea:\n";
            print LIB ".\n";
            print LIB "}\n";
        }

        $dist{$d} = $distUID;
        $distUID++;
    }

    foreach my $xmlfile (@libfiles) {
        my $prefix;

        if ($xmlfile =~ m/xml.(.*\.\d+).bz2$/) {
            $prefix = $1;
        }

        if ($xmlfile =~ m/xml.(.*\.\d+).gz$/) {
            $prefix = $1;
        }

        if ($xmlfile =~ m/xml.(.*\.\d+)$/) {
            $prefix = $1;
        }

        open(L, "> $prefix.frglib.tmp") or die "Failed to open '$prefix.frglib.tmp' for write.\n";
        open(F, "< $prefix.lib")        or die "Failed to open '$prefix.lib' for read\n";

        while (<F>) {
            chomp;
            my ($id, $template, $end, $lib, $mn, $sd, $clr, $clv, $clq) = split '\s+', $_;

            #  Remember which library this frag came from.
            #
            print L "$id\t$dist{$lib}\n";

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

                    #  Uncomment this to do a little sanity checking here
                    #  (uses lots of memory).  It will probably report if
                    #  a template has more than two frags, because the end
                    #  test above will fail.  ????????

                    delete $link{$template}
                } else {
                    $link{$template} = "$id\0$end";
                }
            }

            $numFrag++;
        }
        close(F);
        close(L);

        rename "$prefix.frglib.tmp", "$prefix.frglib";

        my $percMated = 0;
        if ($numFrag > 0) {
            $percMated = int(10000 * 2 * $numMate / $numFrag) / 100;
        }
        print STDERR "$prefix: frags=$numFrag links=$numMate ($percMated%)\n";
    }

    close(LKG);
    close(LIB);

    rename "$prefix.1.lib.frg.tmp", "$prefix.1.lib.frg";
    rename "$prefix.3.lkg.frg.tmp", "$prefix.3.lkg.frg";

    exit(0);
}


sub runFRG ($) {
    my $frgfile = shift @_;
    my $prefix;
    my $index;

    if ($frgfile =~ m/xml.(.*)\.(\d+).bz2$/) {
        $prefix = $1;
        $index  = $2;
        open(X, "bzip2 -dc xml.$prefix.$index.bz2   | $tcat") or die "Failed to open 'xml.$prefix.$index.bz2' for read\n";
        open(F, "bzip2 -dc fasta.$prefix.$index.bz2 | $tcat") or die "Failed to open 'fasta.$prefix.$index.bz2' for read\n";
        open(Q, "bzip2 -dc qual.$prefix.$index.bz2  | $tcat") or die "Failed to open 'qual.$prefix.$index.bz2' for read\n";
        open(L, "< $prefix.$index.frglib")                    or die "Failed to open '$prefix.$index.frglib' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+).gz$/) {
        $prefix = $1;
        $index  = $2;
        open(X, "gzip -dc xml.$prefix.$index.gz   | $tcat") or die "Failed to open 'xml.$prefix.$index.gz' for read\n";
        open(F, "gzip -dc fasta.$prefix.$index.gz | $tcat") or die "Failed to open 'fasta.$prefix.$index.gz' for read\n";
        open(Q, "gzip -dc qual.$prefix.$index.gz  | $tcat") or die "Failed to open 'qual.$prefix.$index.gz' for read\n";
        open(L, "< $prefix.$index.frglib")                  or die "Failed to open '$prefix.$index.frglib' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+)$/) {
        $prefix = $1;
        $index  = $2;
        open(X, "< xml.$prefix.$index")    or die "Failed to open 'xml.$prefix.$index' for read\n";
        open(F, "< fasta.$prefix.$index")  or die "Failed to open 'fasta.$prefix.$index' for read\n";
        open(Q, "< qual.$prefix.$index")   or die "Failed to open 'qual.$prefix.$index' for read\n";
        open(L, "< $prefix.$index.frglib") or die "Failed to open '$prefix.$index.frglib' for read\n";
    }

    if (-e "$prefix.2.$index.frg.bz2") {
        exit(0);
    }

    open(FRG, "| $tcat bzip2 -9c > $prefix.2.$index.frg.bz2.tmp") or die "Failed to open '$prefix.2.$index.frg.bz2.tmp' for write\n";
    #open(FRG, "> $prefix.2.$index.frg.tmp") or die "Failed to open $prefix.2.$index.frg.tmp\n";
    if ($version > 1) {
        print FRG "{VER\n";
        print FRG "ver:$version\n";
        print FRG "}\n";
    }

    my $haveMore = 1;
    $haveMore  = !eof(X);
    $haveMore &= !eof(F);
    $haveMore &= !eof(Q);
    $haveMore &= !eof(L);

    while ($haveMore) {
        my ($xid, $type, $template, $end, $lib, $libsize, $libstddev, $clr, $clv, $clq) = readXML();
        my ($sid, $seq) = readFasta(1);
        my ($qid, $qlt) = readQual(1);

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
                    print FRG "clq:$clq\n" if (defined($clq));
                    print FRG "clr:$clr\n";
                    print FRG "}\n";
                }
            } else {
                print STDERR "ID mismatch: X='$xid' =?= S='$sid'\n";
                print STDERR "ID mismatch: X='$xid' =?= Q='$qid'\n";
                print STDERR "ID mismatch: X='$xid' =?= L='$lid'\n";
                die;
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
    close(FRG);

    rename "$prefix.2.$index.frg.bz2.tmp", "$prefix.2.$index.frg.bz2";

    exit(0);
}



sub runNBL ($) {
    my $frgfile = shift @_;
    my $prefix;
    my $index;

    if ($frgfile =~ m/xml.(.*)\.(\d+).bz2$/) {
        $prefix = $1;
        $index  = $2;
        open(X, "bzip2 -dc xml.$prefix.$index.bz2   | $tcat") or die "Failed to open 'xml.$prefix.$index.bz2' for read\n";
        open(F, "bzip2 -dc fasta.$prefix.$index.bz2 | $tcat") or die "Failed to open 'fasta.$prefix.$index.bz2' for read\n";
        open(Q, "bzip2 -dc qual.$prefix.$index.bz2  | $tcat") or die "Failed to open 'qual.$prefix.$index.bz2' for read\n";
        open(L, "< $prefix.$index.frglib")                    or die "Failed to open '$prefix.$index.frglib' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+).gz$/) {
        $prefix = $1;
        $index  = $2;
        open(X, "gzip -dc xml.$prefix.$index.gz   | $tcat") or die "Failed to open 'xml.$prefix.$index.gz' for read\n";
        open(F, "gzip -dc fasta.$prefix.$index.gz | $tcat") or die "Failed to open 'fasta.$prefix.$index.gz' for read\n";
        open(Q, "gzip -dc qual.$prefix.$index.gz  | $tcat") or die "Failed to open 'qual.$prefix.$index.gz' for read\n";
        open(L, "< $prefix.$index.frglib")                  or die "Failed to open '$prefix.$index.frglib' for read\n";
    }

    if ($frgfile =~ m/xml.(.*)\.(\d+)$/) {
        $prefix = $1;
        $index  = $2;
        open(X, "< xml.$prefix.$index")    or die "Failed to open 'xml.$prefix.$index' for read\n";
        open(F, "< fasta.$prefix.$index")  or die "Failed to open 'fasta.$prefix.$index' for read\n";
        open(Q, "< qual.$prefix.$index")   or die "Failed to open 'qual.$prefix.$index' for read\n";
        open(L, "< $prefix.$index.frglib") or die "Failed to open '$prefix.$index.frglib' for read\n";
    }

    if (-e "$prefix.2.$index.frg.bz2") {
        exit(0);
    }

    open(FNA, "> $prefix.2.$index.fna.tmp")      or die "Failed to open '$prefix.2.$index.fna.tmp' for write\n";
    open(QNA, "> $prefix.2.$index.fna.qual.tmp") or die "Failed to open '$prefix.2.$index.fna.qual.tmp' for write\n";

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

    rename "$prefix.2.$index.fna.tmp",      "$prefix.2.$index.fna";
    rename "$prefix.2.$index.fna.qual.tmp", "$prefix.2.$index.fna.qual";

    exit(0);
}
