
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
 #  This file is derived from:
 #
 #    src/pipelines/ca3g/Defaults.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-21
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-19
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Defaults;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getCommandLineOptions addCommandLineOption addCommandLineError writeLog getNumberOfCPUs getPhysicalMemorySize diskSpace printOptions printHelp printCitation addSequenceFile setParametersFromFile setParametersFromCommandLine checkJava checkGnuplot checkParameters getGlobal setGlobal setGlobalIfUndef setDefaults setVersion);

use strict;
use Cwd qw(getcwd abs_path);
use Carp qw(cluck);
use Sys::Hostname;
use Text::Wrap;
use File::Basename;   #  dirname

my %global;    #  Parameter value
my %synops;    #  Parameter description (for -defaults)
my %synnam;    #  Parameter name (beacuse the key is lowercase)

my $cLineOpts = undef;
my $specLog   = "";



sub getGlobal ($) {
    my $var = shift @_;

    $var =~ tr/A-Z/a-z/;

    #  We lost the use of caFailure in Defaults.pm (because it was moved to
    #  Execution.pm so it can run stuff) here, so duplicate the functionality.
    #  This should only trigger on static pipeline errors (i.e., no depending
    #  on reads input) and so should never occur in the wild.

    if (!exists($global{$var})) {
        print STDERR "================================================================================\n";
        print STDERR "Unknown parameter '$var' accessed.  Stack trace:\n";
        cluck;
        exit(1);
    }

    return($global{$var});
}


sub setGlobalSpecialization ($@) {
    my $val = shift @_;

    foreach my $var (@_) {
        $global{$var} = $val;
    }

    return(1);
}


sub setGlobal ($$) {
    my $VAR = shift @_;
    my $var = $VAR;
    my $val = shift @_;
    my $set = 0;

    $var =~ tr/A-Z/a-z/;

    $val = undef  if ($val eq "undef");   #  Set to undefined, the default for many of the options.
    $val = undef  if ($val eq "");

    #  Map 'true'/'false' et al. to 0/1.

    $val = 0  if (($val =~ m/^false$/i) || ($val =~ m/^f$/i));
    $val = 1  if (($val =~ m/^true$/i)  || ($val =~ m/^t$/i));

    #  Grid options

    foreach my $opt ("gridoptions") {
        $set += setGlobalSpecialization($val, ("${opt}corovl",  "${opt}obtovl",  "${opt}utgovl"))   if ($var eq "${opt}ovl");
        $set += setGlobalSpecialization($val, ("${opt}cormhap", "${opt}obtmhap", "${opt}utgmhap"))  if ($var eq "${opt}mhap");
        $set += setGlobalSpecialization($val, ("${opt}cormmap", "${opt}obtmmap", "${opt}utgmmap"))  if ($var eq "${opt}mmap");
    }

    foreach my $opt ("memory",
                     "threads",
                     "concurrency") {
        $set += setGlobalSpecialization($val, ( "corovl${opt}",  "obtovl${opt}",  "utgovl${opt}"))  if ($var eq  "ovl${opt}");
        $set += setGlobalSpecialization($val, ("cormhap${opt}", "obtmhap${opt}", "utgmhap${opt}"))  if ($var eq "mhap${opt}");
        $set += setGlobalSpecialization($val, ("cormmap${opt}", "obtmmap${opt}", "utgmmap${opt}"))  if ($var eq "mmap${opt}");
    }

    #  Overlapping algorithm choice options

    foreach my $opt ("overlapper",
                     "realign") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  OverlapInCore options

    foreach my $opt ("ovlerrorrate",
                     "ovlhashblocklength",
                     "ovlrefblocksize",
                     "ovlrefblocklength",
                     "ovlhashbits",
                     "ovlhashload",
                     "ovlmersize",
                     "ovlmerthreshold",
                     "ovlmerdistinct",
                     "ovlmertotal",
                     "ovlfrequentmers") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  Mhap options

    foreach my $opt ("mhapblocksize",
                     "mhapmersize",
                     "mhapsensitivity",
                     "mhapfilterunique",
                     "mhapfilterthreshold",
                     "mhapnotf") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  MiniMap options

    foreach my $opt ("mmapblocksize",
                     "mmapmersize") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  Handle the two error rate aliases.

    if ($var eq "rawerrorrate") {
        setGlobalIfUndef("corOvlErrorRate", $val);
        setGlobalIfUndef("corErrorRate",    $val);
        return;
    }

    if ($var eq "correctederrorrate") {
        setGlobalIfUndef("obtOvlErrorRate", $val);
        setGlobalIfUndef("obtErrorRate",    $val);
        setGlobalIfUndef("utgOvlErrorRate", $val);
        setGlobalIfUndef("utgErrorRate",    $val);
        setGlobalIfUndef("cnsErrorRate",    $val);
        return;
    }

    return  if ($set > 0);

    #if ($var eq "canuiteration") {
    #    print STDERR "-- WARNING: set canuIteration to $val\n";
    #}

    #  If we get a parameter we don't understand, we should be parsing command line options or
    #  reading spec files, and we can let the usual error handling handle it.

    addCommandLineError("ERROR:  Parameter '$VAR' is not known.\n")   if (!exists($global{$var}));

    $global{$var} = $val;
}



sub setGlobalIfUndef ($$) {
    my $var = shift @_;
    my $val = shift @_;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    return  if (defined($global{$var}));

    $global{$var} = $val;
}



sub getCommandLineOptions () {
    return($cLineOpts);
}



sub addCommandLineOption ($) {
    my $opt = shift @_;

    return   if ($opt =~ m/canuIteration=/);   #  Ignore canu resetting canuIteration

    $cLineOpts .= " "   if (defined($cLineOpts) && ($cLineOpts !~ m/\s$/));
    $cLineOpts .= $opt;
}



sub addCommandLineError($) {
    $global{'errors'} .= shift @_;
}



sub writeLog () {
    my $time = time();
    my $host = hostname();
    my $pid  = $$;

    open(F, "> canu-logs/${time}_${host}_${pid}_canu");
    print F $specLog;
    close(F);
}



#
#  Host management - these really belong in 'Execution.pm' (or 'Utilities.pm') but can't go there
#  (Execution.pm) and be used here too.
#

sub getNumberOfCPUs () {
    my $os   = $^O;
    my $ncpu = 1;

    #  See http://stackoverflow.com/questions/6481005/obtain-the-number-of-cpus-cores-in-linux

    if ($os eq "freebsd") {
        $ncpu = int(`/sbin/sysctl -n hw.ncpu`);
    }

    if ($os eq "darwin") {
        $ncpu = int(`/usr/bin/getconf _NPROCESSORS_ONLN`);
    }

    if ($os eq "linux" || $os eq "cygwin") {
        $ncpu = int(`getconf _NPROCESSORS_ONLN`);
    }

    return($ncpu);
}


sub getPhysicalMemorySize () {
    my $os     = $^O;
    my $memory = 1;

    if ($os eq "freebsd") {
        $memory = `/sbin/sysctl -n hw.physmem` / 1024 / 1024 / 1024;
    }

    if ($os eq "darwin") {
        $memory = `/usr/sbin/sysctl -n hw.memsize` / 1024 / 1024 / 1024;
    }

    if ($os eq "linux" || $os eq "cygwin") {
        open(F, "< /proc/meminfo");        #  Way to go, Linux!  Make it easy on us!
        while (<F>) {
            if (m/MemTotal:\s+(\d+)/) {
                $memory = $1 / 1024 / 1024;
            }
        }
        close(F);
    }

    return(int($memory + 0.5));  #  Poor man's rounding
}



sub diskSpace ($) {
    my  $dir                          = dirname($_[0]);
    my ($total, $used, $free, $avail) = (0, 0, 0, 0);

    if (-d $dir) {
        open(DF, "df -P -k $dir |");
        while (<DF>) {
            chomp;

            if (m/^(.*)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+%)\s+(.*)$/) {
                $total = int($2 / 1048.576) / 1000;
                $used  = int($3 / 1048.576) / 1000;
                $free  = int($4 / 1048.576) / 1000;
                $avail = int($4 / 1048.576) / 1000;  #  Possibly limited by quota?
            }
        }
        close(DF);
    }

    #print STDERR "Disk space: total $total GB, used $used GB, free $free GB, available $avail GB\n";

    return (wantarray) ? ($total, $used, $free, $avail) : $avail;
}


sub printOptions () {
    my $pretty = 0;

    foreach my $k (sort values %synnam) {
        my $o = $k;
        my $u = $synops{$k};

        next   if (length($u) == 0);

        if ($pretty == 0) {
            $o = substr("$k                                    ", 0, 40);

        } else {
            $Text::Wrap::columns = 60;

            $o = "$o\n";
            $u = wrap("    ", "    ", $u) . "\n";
        }

        print "$o$u\n";
    }
}


sub printHelp (@) {
    my $force = shift @_;

    return   if (!defined($force) && !defined($global{"errors"}));

    print "\n";
    print "usage:   canu [-version] [-citation] \\\n";
    print "              [-correct | -trim | -assemble | -trim-assemble] \\\n";
    print "              [-s <assembly-specifications-file>] \\\n";
    print "               -p <assembly-prefix> \\\n";
    print "               -d <assembly-directory> \\\n";
    print "               genomeSize=<number>[g|m|k] \\\n";
    print "              [other-options] \\\n";
    print "              [-pacbio-raw |\n";
    print "               -pacbio-corrected |\n";
    print "               -nanopore-raw |\n";
    print "               -nanopore-corrected] file1 file2 ...\n";
    print "\n";
    print "example: canu -d run1 -p godzilla genomeSize=1g -nanopore-raw reads/*.fasta.gz \n";
    print "\n";
    print "\n";
    print "  To restrict canu to only a specific stage, use:\n";
    print "    -correct       - generate corrected reads\n";
    print "    -trim          - generate trimmed reads\n";
    print "    -assemble      - generate an assembly\n";
    print "    -trim-assemble - generate trimmed reads and then assemble them\n";
    print "\n";
    print "  The assembly is computed in the -d <assembly-directory>, with output files named\n";
    print "  using the -p <assembly-prefix>.  This directory is created if needed.  It is not\n";
    print "  possible to run multiple assemblies in the same directory.\n";
    print "\n";
    print "  The genome size should be your best guess of the haploid genome size of what is being\n";
    print "  assembled.  It is used primarily to estimate coverage in reads, NOT as the desired\n";
    print "  assembly size.  Fractional values are allowed: '4.7m' equals '4700k' equals '4700000'\n";
    print "\n";
    print "  Some common options:\n";
    print "    useGrid=string\n";
    print "      - Run under grid control (true), locally (false), or set up for grid control\n";
    print "        but don't submit any jobs (remote)\n";
    print "    rawErrorRate=fraction-error\n";
    print "      - The allowed difference in an overlap between two raw uncorrected reads.  For lower\n";
    print "        quality reads, use a higher number.  The defaults are 0.300 for PacBio reads and\n";
    print "        0.500 for Nanopore reads.\n";
    print "    correctedErrorRate=fraction-error\n";
    print "      - The allowed difference in an overlap between two corrected reads.  Assemblies of\n";
    print "        low coverage or data with biological differences will benefit from a slight increase\n";
    print "        in this.  Defaults are 0.045 for PacBio reads and 0.144 for Nanopore reads.\n";
    print "    gridOptions=string\n";
    print "      - Pass string to the command used to submit jobs to the grid.  Can be used to set\n";
    print "        maximum run time limits.  Should NOT be used to set memory limits; Canu will do\n";
    print "        that for you.\n";
    print "    minReadLength=number\n";
    print "      - Ignore reads shorter than 'number' bases long.  Default: 1000.\n";
    print "    minOverlapLength=number\n";
    print "      - Ignore read-to-read overlaps shorter than 'number' bases long.  Default: 500.\n";
    print "  A full list of options can be printed with '-options'.  All options can be supplied in\n";
    print "  an optional sepc file with the -s option.\n";
    print "\n";
    print "  Reads can be either FASTA or FASTQ format, uncompressed, or compressed with gz, bz2 or xz.\n";
    print "  Reads are specified by the technology they were generated with:\n";
    print "    -pacbio-raw         <files>\n";
    print "    -pacbio-corrected   <files>\n";
    print "    -nanopore-raw       <files>\n";
    print "    -nanopore-corrected <files>\n";
    print "\n";
    print "Complete documentation at http://canu.readthedocs.org/en/latest/\n";
    print "\n";

    if (defined($global{'errors'})) {
        print "$global{'errors'}";
        print "\n";
    }

    exit(1);
}


sub printCitation ($) {
    my $prefix = shift @_;

    print STDERR "${prefix}Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM.\n";
    print STDERR "${prefix}Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation.\n";
    print STDERR "${prefix}Genome Res. 2017 May;27(5):722-736.\n";
    print STDERR "${prefix}http://doi.org/10.1101/gr.215087.116\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}Read and contig alignments during correction, consensus and GFA building use:\n";
    print STDERR "${prefix}  Šošic M, Šikic M.\n";
    print STDERR "${prefix}  Edlib: a C/C ++ library for fast, exact sequence alignment using edit distance.\n";
    print STDERR "${prefix}  Bioinformatics. 2017 May 1;33(9):1394-1395.\n";
    print STDERR "${prefix}  http://doi.org/10.1093/bioinformatics/btw753\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}Overlaps are generated using:\n";
    print STDERR "${prefix}  Berlin K, et al.\n";
    print STDERR "${prefix}  Assembling large genomes with single-molecule sequencing and locality-sensitive hashing.\n";
    print STDERR "${prefix}  Nat Biotechnol. 2015 Jun;33(6):623-30.\n";
    print STDERR "${prefix}  http://doi.org/10.1038/nbt.3238\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}  Myers EW, et al.\n";
    print STDERR "${prefix}  A Whole-Genome Assembly of Drosophila.\n";
    print STDERR "${prefix}  Science. 2000 Mar 24;287(5461):2196-204.\n";
    print STDERR "${prefix}  http://doi.org/10.1126/science.287.5461.2196\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}  Li H.\n";
    print STDERR "${prefix}  Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences.\n";
    print STDERR "${prefix}  Bioinformatics. 2016 Jul 15;32(14):2103-10.\n";
    print STDERR "${prefix}  http://doi.org/10.1093/bioinformatics/btw152\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}Corrected read consensus sequences are generated using an algorithm derived from FALCON-sense:\n";
    print STDERR "${prefix}  Chin CS, et al.\n";
    print STDERR "${prefix}  Phased diploid genome assembly with single-molecule real-time sequencing.\n";
    print STDERR "${prefix}  Nat Methods. 2016 Dec;13(12):1050-1054.\n";
    print STDERR "${prefix}  http://doi.org/10.1038/nmeth.4035\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}Contig consensus sequences are generated using an algorithm derived from pbdagcon:\n";
    print STDERR "${prefix}  Chin CS, et al.\n";
    print STDERR "${prefix}  Nonhybrid, finished microbial genome assemblies from long-read SMRT sequencing data.\n";
    print STDERR "${prefix}  Nat Methods. 2013 Jun;10(6):563-9\n";
    print STDERR "${prefix}  http://doi.org/10.1038/nmeth.2474\n";
    print STDERR "${prefix}\n";
}




sub makeAbsolute ($) {
    my $var = shift @_;
    my $val = getGlobal($var);
    my $abs = abs_path($val);

    if (defined($val) && ($val != $abs)) {
        setGlobal($var, $abs);
        $val =~ s/\\\"/\"/g;
        $val =~ s/\"/\\\"/g;
        $val =~ s/\\\$/\$/g;
        $val =~ s/\$/\\\$/g;

        addCommandLineOption("'$var=$val'");
    }
}



sub fixCase ($) {
    my $var = shift @_;
    my $val = getGlobal($var);

    if (defined($val)) {
        $val =~ tr/A-Z/a-z/;
        setGlobal($var, $val);
    }
}



sub addSequenceFile ($$@) {
    my $dir   = shift @_;
    my $file  = shift @_;
    my $err   = shift @_;

    return(undef)             if (!defined($file));   #  No file name?  Nothing to do.
    $file = "$dir/$file"      if (defined($dir));     #  If $dir defined, assume file is in there.
    return($file)             if ($file =~ m!/!);     #  If already a full path, use that.
    return(abs_path($file))   if (-e $file);          #  If found, return the full path.

    #  And if not found, report an error, unless told not to.  This is because on the command
    #  line, the first word after -pacbio-raw must exist, but all the other words could
    #  be files or options.

    addCommandLineError("ERROR: Input read file '$file' not found.\n")  if (defined($err));

    return(undef);
}



sub setParametersFromFile ($$@) {
    my $specFile  = shift @_;
    my $readdir   = shift @_;
    my @fragFiles = @_;

    #  Client should be ensuring that the file exists before calling this function.
    die "specFile '$specFile' not found.\n"  if (! -e "$specFile");

    $specLog .= "\n";
    $specLog .= "###\n";
    $specLog .= "###  Reading options from '$specFile'\n";
    $specLog .= "###\n";
    $specLog .= "\n";

    #  We lost the use of caExit() here (moved to Execution.pm) and so can't call it.
    #  Just die.

    open(F, "< $specFile") or die("can't open '$specFile' for reading: $!\n");

    while (<F>) {
        $specLog .= $_;

        s/^\s+//;
        s/\s+$//;

        next if (m/^#/);
        next if (length($_) eq 0);

        #  First, figure out the two words.

        my $one;
        my $two;
        my $opt;

        if (m/^-(pacbio|nanopore)-(corrected|raw)\s+(.*)\s*$/) {   #  Comments not allowed, because then we can't decide
            $one  = "-$1-$2";                                      #  if the # is a comment, or part of the file!
            $two = $3;                                             #  e.g.,   this_is_file_#1   vs
            $opt = 0;                                              #          this_is_the_only_file#no more data
        }

        elsif (m/^(\w*)\s*=\s*([^#]*)\s*#*.*?$/) {   #  Word two won't match a #, but will gobble up spaces at the end.
            $one = $1;                               #  Then, we can match a #, and any amount of comment, minimally.
            $two = $2;                               #  If word two is made non-greedy, it will shrink to nothing, as
            $opt = 1;                                #  the last bit will gobble up everything, since we're allowed
        }                                            #  to match zero #'s in between.

        else {
            addCommandLineError("ERROR:  File not found or unknown specFile option line '$_'.\n");
        }

        #  Now, clean up the second word to handle quotes.

        $two =~ s/^\s+//;   #  There can be spaces from the greedy match.
        $two =~ s/\s+$//;

        $two = $1   if ($two =~ m/^'(.+)'$/);    #  Remove single quotes   |   But don't allowed mixed quotes; users
        $two = $1   if ($two =~ m/^"(.+)"$/);    #  Remove double quotes   |   should certainly know better

        #  And do something.

        if ($opt == 1) {
            $two =~ s/^\s+//;  #  Remove spaces again.  They'll just confuse our option processing.
            $two =~ s/\s+$//;

            setGlobal($one, $two);
        }

        else {
            my $file = addSequenceFile($readdir, $two, 1);   #  Don't remove spaces.  File could be " file ", for some stupid reason.

            if (defined($file)) {
                push @fragFiles, "$one\0$file";
            } else {
                addCommandLineError("ERROR:  File not found in spec file option '$_'\n");
            }
        }
    }
    close(F);

    return(@fragFiles);
}



sub setParametersFromCommandLine(@) {
    my @specOpts = @_;

    if (scalar(@specOpts) > 0) {
        $specLog .= "\n";
        $specLog .= "###\n";
        $specLog .= "###  Reading options from the command line.\n";
        $specLog .= "###\n";
        $specLog .= "\n";
    }

    foreach my $s (@specOpts) {
        $specLog .= "$s\n";

        if ($s =~ m/\s*(\w*)\s*=(.*)/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            setGlobal($var, $val);
        } else {
            addCommandLineError("ERROR:  Misformed command line option '$s'.\n");
        }
    }
}



sub setExecDefaults ($$) {
    my $tag         = shift @_;
    my $name        = shift @_;

    $global{"gridOptions${tag}"}   = undef;
    $synops{"gridOptions${tag}"}   = "Grid engine options applied to $name jobs";

    $global{"${tag}Memory"}        = undef;
    $synops{"${tag}Memory"}        = "Amount of memory, in gigabytes, to use for $name jobs";

    $global{"${tag}Threads"}       = undef;
    $synops{"${tag}Threads"}       = "Number of threads to use for $name jobs";

    $global{"${tag}StageSpace"}    = undef;
    $synops{"${tag}StageSpace"}    = "Amount of local disk space needed to stage data for $name jobs";

    $global{"${tag}Concurrency"}   = undef;
    $synops{"${tag}Concurrency"}   = "If grid not enabled, number of $name jobs to run at the same time; default is n_proc / n_threads";
}



sub setOverlapDefault ($$$$) {
    my $tag         = shift @_;
    my $var         = shift @_;
    my $value       = shift @_;
    my $description = shift @_;

    $global{"${tag}${var}"}  = $value;
    $synops{"${tag}${var}"}  = $description;
    $synops{      "${var}"}  = $description;
}



sub setOverlapDefaults ($$$) {
    my $tag     = shift @_;  #  If 'cor', some parameters are loosened for raw pacbio reads
    my $name    = shift @_;
    my $default = shift @_;  #  Sets ${tag}Overlapper

    #  Which overlapper to use.

    setOverlapDefault($tag, "Overlapper",          $default,                  "Which overlap algorithm to use for $name");
    setOverlapDefault($tag, "ReAlign",             0,                         "Refine overlaps by computing the actual alignment: 'true' or 'false'.  Not useful for overlapper=ovl.  Uses ${tag}OvlErrorRate");

    #  OverlapInCore parameters.

    setOverlapDefault($tag, "OvlHashBlockLength",  undef,                     "Amount of sequence (bp) to load into the overlap hash table");
    setOverlapDefault($tag, "OvlRefBlockSize",     undef,                     "Number of reads to search against the hash table per batch");
    setOverlapDefault($tag, "OvlRefBlockLength",   0,                         "Amount of sequence (bp) to search against the hash table per batch");
    setOverlapDefault($tag, "OvlHashBits",         ($tag eq "cor") ? 18 : 23, "Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per ${tag}OvlHashBlockLength");
    setOverlapDefault($tag, "OvlHashLoad",         0.75,                      "Maximum hash table load.  If set too high, table lookups are inefficent; if too low, search overhead dominates run time; default 0.75");
    setOverlapDefault($tag, "OvlMerSize",          ($tag eq "cor") ? 19 : 22, "K-mer size for seeds in overlaps");
    setOverlapDefault($tag, "OvlMerThreshold",     "auto",                    "K-mer frequency threshold; mers more frequent than this count are ignored; default 'auto'");
    setOverlapDefault($tag, "OvlMerDistinct",      undef,                     "K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps");
    setOverlapDefault($tag, "OvlMerTotal",         undef,                     "K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps");
    setOverlapDefault($tag, "OvlFrequentMers",     undef,                     "Do not seed overlaps with these kmers (fasta format)");
    setOverlapDefault($tag, "OvlFilter",           undef,                     "Filter overlaps based on expected kmers vs observed kmers");

    #  Mhap parameters.  FilterThreshold MUST be a string, otherwise it gets printed in scientific notation (5e-06) which java doesn't understand.

    setOverlapDefault($tag, "MhapVersion",         "2.1.2",                   "Version of the MHAP jar file to use");
    setOverlapDefault($tag, "MhapFilterThreshold", "0.000005",                "Value between 0 and 1. kmers which comprise more than this percentage of the input are downweighted");
    setOverlapDefault($tag, "MhapFilterUnique",    undef,                     "Expert option: True or false, supress the low-frequency k-mer distribution based on them being likely noise and not true overlaps. Threshold auto-computed based on error rate and coverage.");
    setOverlapDefault($tag, "MhapNoTf",            undef,                     "Expert option: True or false, do not use tf weighting, only idf of tf-idf.");
    setOverlapDefault($tag, "MhapOptions",         undef,                     "Expert option: free-form parameters to pass to MHAP.");
    setOverlapDefault($tag, "MhapBlockSize",       3000,                      "Number of reads per 1GB; memory * blockSize = the size of  block loaded into memory per job");
    setOverlapDefault($tag, "MhapMerSize",         ($tag eq "cor") ? 16 : 16, "K-mer size for seeds in mhap");
    setOverlapDefault($tag, "MhapOrderedMerSize",  ($tag eq "cor") ? 12 : 18, "K-mer size for second-stage filter in mhap");
    setOverlapDefault($tag, "MhapSensitivity",     undef,                     "Coarse sensitivity level: 'low', 'normal' or 'high'.  Set automatically based on coverage; 'high' <= 30x < 'normal' < 60x <= 'low'");

    #  MiniMap parameters.

    setOverlapDefault($tag, "MMapBlockSize",       6000,                      "Number of reads per 1GB; memory * blockSize = the size of  block loaded into memory per job");
    setOverlapDefault($tag, "MMapMerSize",         ($tag eq "cor") ? 15 : 21, "K-mer size for seeds in minmap");
}



sub setDefault ($$$) {
    my $var         = shift @_;
    my $value       = shift @_;
    my $description = shift @_;

    $global{$var} = $value;
    $synops{$var} = $description;
}



sub setDefaults () {

    #####  Internal stuff

    $global{"errors"}                      = undef;   #  Command line errors
    $global{"version"}                     = undef;   #  Reset at the end of this function, once we know where binaries are.
    $global{"availablehosts"}              = undef;   #  Internal list of cpus-memory-nodes describing the grid

    $global{"canuiteration"}               = 0;
    $global{"canuiterationmax"}            = 2;

    $global{"onexitdir"}                   = undef;   #  Copy of $wrk, for caExit() and caFailure() ONLY.
    $global{"onexitnam"}                   = undef;   #  Copy of $asm, for caExit() and caFailure() ONLY.

    #####  Meta options (no $global for these, only synopsis), more of these, many many more, are defined in setOverlapDefaults().

    $synops{"rawErrorRate"}                = "Expected fraction error in an alignment of two uncorrected reads";
    $synops{"correctedErrorRate"}          = "Expected fraction error in an alignment of two corrected reads";

    #####  General Configuration Options (aka miscellany)

    my $java = (exists $ENV{"JAVA_HOME"} && -e "$ENV{'JAVA_HOME'}/bin/java") ? "$ENV{'JAVA_HOME'}/bin/java" : "java";

    setDefault("showNext",            undef,     "Don't run any commands, just report what would run");
    setDefault("pathMap",             undef,     "File with a hostname to binary directory map; binary directories must be absolute paths");
    setDefault("shell",               "/bin/sh", "Command interpreter to use; sh-compatible (e.g., bash), NOT C-shell (csh or tcsh); default '/bin/sh'");
    setDefault("java",                $java,     "Java interpreter to use; at least version 1.8; default 'java'");
    setDefault("gnuplot",             "gnuplot", "Path to the gnuplot executable");
    setDefault("gnuplotImageFormat",  undef,     "Image format that gnuplot will generate, used in HTML reports.  Default: based on gnuplot, 'png', 'svg' or 'gif'");
    setDefault("gnuplotTested",       0,         "If set, skip the initial testing of gnuplot");
    setDefault("stageDirectory",      undef,     "If set, copy heavily used data to this node-local location");

    #####  Cleanup and Termination options

    setDefault("saveOverlaps",        0,     "Save intermediate overlap files, almost never a good idea");
    setDefault("saveReadCorrections", 0,     "Save intermediate read correction files, almost never a good idea");
    setDefault("saveMerCounts",       0,     "Save full mer counting results, sometimes useful");
    setDefault("onSuccess",           undef, "Full path to command to run on successful completion");
    setDefault("onFailure",           undef, "Full path to command to run on failure");

    #####  Error Rates

    setDefault("corOvlErrorRate",     undef, "Overlaps above this error rate are not computed");
    setDefault("obtOvlErrorRate",     undef, "Overlaps at or below this error rate are used to trim reads");
    setDefault("utgOvlErrorRate",     undef, "Overlaps at or below this error rate are used to trim reads");
    setDefault("utgErrorRate",        undef, "Overlaps at or below this error rate are used to construct contigs");
    setDefault("utgGraphDeviation",   6,     "Overlaps this much above median will not be used for initial graph construction");
    setDefault("utgRepeatDeviation",  3,     "Overlaps this much above mean unitig error rate will not be used for repeat splitting");
    setDefault("utgRepeatConfusedBP", 2100,  "Repeats where the next best edge is at least this many bp shorter will not be split");
    setDefault("corErrorRate",        undef, "Only use raw alignments below this error rate to construct corrected reads");
    setDefault("cnsErrorRate",        undef, "Consensus expects alignments at about this error rate");

    #####  Minimums and maximums

    setDefault("minReadLength",    1000, "Reads shorter than this length are not loaded into the assembler; default 1000");
    setDefault("minOverlapLength",  500, "Overlaps shorter than this length are not computed; default 500");

    setDefault("minMemory",        undef, "Minimum amount of memory needed to compute the assembly (do not set unless prompted!)");
    setDefault("maxMemory",        undef, "Maximum memory to use by any component of the assembler");

    setDefault("minThreads",       undef, "Minimum number of compute threads suggested to compute the assembly");
    setDefault("maxThreads",       undef, "Maximum number of compute threads to use by any component of the assembler");

    #####  Stopping conditions

    setDefault("stopOnReadQuality", 1,     "Stop if a significant portion of the input data is too short or has quality value or base composition errors");
    setDefault("stopAfter",         undef, "Stop after a specific algorithm step is completed");

    #####  Grid Engine configuration, internal parameters.  These are filled out in canu.pl, right after this function returns.

    setDefault("gridEngine",                          undef, "Grid engine configuration, not documented");
    setDefault("gridEngineSubmitCommand",             undef, "Grid engine configuration, not documented");
    setDefault("gridEngineNameOption",                undef, "Grid engine configuration, not documented");
    setDefault("gridEngineArrayOption",               undef, "Grid engine configuration, not documented");
    setDefault("gridEngineArrayName",                 undef, "Grid engine configuration, not documented");
    setDefault("gridEngineArrayMaxJobs",              undef, "Grid engine configuration, not documented");
    setDefault("gridEngineOutputOption",              undef, "Grid engine configuration, not documented");
    setDefault("gridEnginePropagateCommand",          undef, "Grid engine configuration, not documented");
    setDefault("gridEngineThreadsOption",             undef, "Grid engine configuration, not documented");
    setDefault("gridEngineMemoryOption",              undef, "Grid engine configuration, not documented");
    setDefault("gridEngineMemoryUnits",               undef, "Grid engine configuration, not documented");
    setDefault("gridEngineNameToJobIDCommand",        undef, "Grid engine configuration, not documented");
    setDefault("gridEngineNameToJobIDCommandNoArray", undef, "Grid engine configuration, not documented");
    setDefault("gridEngineStageOption",               undef, "Grid engine configuration, not documented");
    setDefault("gridEngineTaskID",                    undef, "Grid engine configuration, not documented");
    setDefault("gridEngineArraySubmitID",             undef, "Grid engine configuration, not documented");
    setDefault("gridEngineJobID",                     undef, "Grid engine configuration, not documented");

    #####  Grid Engine Pipeline

    setDefault("useGrid", 1, "If 'true', enable grid-based execution; if 'false', run all jobs on the local machine; if 'remote', create jobs for grid execution but do not submit; default 'true'");

    foreach my $c (qw(BAT GFA CNS COR MERYL CORMHAP CORMMAP COROVL OBTMHAP OBTMMAP OBTOVL OEA OVB OVS RED UTGMHAP UTGMMAP UTGOVL)) {
        setDefault("useGrid$c", 1, "If 'true', run module $c under grid control; if 'false' run locally.");
    }

    #####  Grid Engine configuration, for each step of the pipeline

    setDefault("gridOptions",           undef,  "Grid engine options applied to all jobs");
    setDefault("gridOptionsExecutive",  undef,  "Grid engine options applied to the canu executive script");
    setDefault("gridOptionsJobName",    undef,  "Grid jobs job-name suffix");

    #####  Grid Engine configuration and parameters, for each step of the pipeline (memory, threads)

    setExecDefaults("meryl",   "mer counting");
    setExecDefaults("cor",     "read correction");

    setExecDefaults("corovl",  "overlaps for correction");
    setExecDefaults("obtovl",  "overlaps for trimming");
    setExecDefaults("utgovl",  "overlaps for unitig construction");

    setExecDefaults("cormhap", "mhap overlaps for correction");
    setExecDefaults("obtmhap", "mhap overlaps for trimming");
    setExecDefaults("utgmhap", "mhap overlaps for unitig construction");

    setExecDefaults("cormmap", "mmap overlaps for correction");
    setExecDefaults("obtmmap", "mmap overlaps for trimming");
    setExecDefaults("utgmmap", "mmap overlaps for unitig construction");

    setExecDefaults("ovb",     "overlap store bucketizing");
    setExecDefaults("ovs",     "overlap store sorting");

    setExecDefaults("red",     "read error detection");
    setExecDefaults("oea",     "overlap error adjustment");

    setExecDefaults("bat",     "unitig construction");
    setExecDefaults("cns",     "unitig consensus");
    setExecDefaults("gfa",     "graph alignment and processing");

    #####  Object Storage

    setDefault("objectStore",          undef,  "Type of object storage used; not ready for production yet");
    setDefault("objectStoreClient",    undef,  "Path to the command line client used to access the object storage");
    setDefault("objectStoreNameSpace", undef,  "Object store parameters; specific to the type of objectStore used");

    #####  Overlapper

    setOverlapDefaults("cor", "correction",             "mhap");  #  Overlaps computed for correction
    setOverlapDefaults("obt", "overlap based trimming", "ovl");   #  Overlaps computed for trimming
    setOverlapDefaults("utg", "unitig construction",    "ovl");   #  Overlaps computed for unitigging

    ##### Overlap Store

    setDefault("ovsMethod", undef, "Use the 'sequential' or 'parallel' algorithm for constructing an overlap store; default 'sequential'");

    #####  Mers

    setDefault("merylMemory",      undef,  "Amount of memory, in gigabytes, to use for mer counting");
    setDefault("merylThreads",     undef,  "Number of threads to use for mer counting");
    setDefault("merylConcurrency", undef,  "Unused, there is only one process");

    #####  Overlap Based Trimming

    setDefault("obtErrorRate",       undef, "Stringency of overlaps to use for trimming");
    setDefault("trimReadsOverlap",   1,     "Minimum overlap between evidence to make contiguous trim; default '1'");
    setDefault("trimReadsCoverage",  1,     "Minimum depth of evidence to retain bases; default '1'");

    #$global{"splitReads..."}               = 1;
    #$synops{"splitReads..."}               = "";

    #####  Fragment/Overlap Error Correction

    setDefault("enableOEA",      1,     "Do overlap error adjustment - comprises two steps: read error detection (RED) and overlap error adjustment (OEA); default 'true'");
    setDefault("redBatchSize",   undef, "Number of reads per fragment error detection batch");
    setDefault("redBatchLength", undef, "Number of bases per fragment error detection batch");
    setDefault("oeaBatchSize",   undef, "Number of reads per overlap error correction batch");
    setDefault("oeaBatchLength", undef, "Number of bases per overlap error correction batch");

    #####  Unitigger & BOG & bogart Options

    setDefault("unitigger",      "bogart", "Which unitig algorithm to use; only 'bogart' supported; default 'bogart'");
    setDefault("genomeSize",     undef, "An estimate of the size of the genome");
    setDefault("batOptions",     undef, "Advanced options to bogart");
    setDefault("batMemory",      undef, "Approximate maximum memory usage, in gigabytes, default is the maxMemory limit");
    setDefault("batThreads",     undef, "Number of threads to use; default is the maxThreads limit");
    setDefault("batConcurrency", undef, "Unused, only one process supported");

    setDefault("contigFilter",   "2 0 1.0 0.5 5",   "Parameters to filter out 'unassembled' unitigs.  Five values: minReads minLength singleReadSpan lowCovFraction lowCovDepth");

    #####  Consensus Options

    setDefault("cnsPartitions",   undef,       "Partition consensus into N jobs");
    setDefault("cnsPartitionMin", undef,       "Don't make a consensus partition with fewer than N reads");
    setDefault("cnsMaxCoverage",  40,          "Limit unitig consensus to at most this coverage; default '0' = unlimited");
    setDefault("cnsConsensus",    "pbdagcon",  "Which consensus algorithm to use; 'pbdagcon' (fast, reliable); 'utgcns' (multialignment output); 'quick' (single read mosaic); default 'pbdagcon'");

    #####  Correction Options

    setDefault("corPartitions",                undef,        "Partition read correction into N jobs");
    setDefault("corPartitionMin",              undef,        "Don't make a read correction partition with fewer than N reads");
    setDefault("corMinEvidenceLength",         undef,        "Limit read correction to only overlaps longer than this; default: unlimited");
    setDefault("corMaxEvidenceErate",          undef,        "Limit read correction to only overlaps at or below this fraction error; default: unlimited");
    setDefault("corMaxEvidenceCoverageGlobal", "1.0x",       "Limit reads used for correction to supporting at most this coverage; default: '1.0x' = 1.0 * estimated coverage");
    setDefault("corMaxEvidenceCoverageLocal",  "2.0x",       "Limit reads being corrected to at most this much evidence coverage; default: '2.0x' = 2.0 * estimated coverage");
    setDefault("corOutCoverage",               40,           "Only correct the longest reads up to this coverage; default 40");
    setDefault("corMinCoverage",               undef,        "Minimum number of bases supporting each corrected base, if less than this sequences are split; default based on input read coverage: 0 <= 30x < 4 < 60x <= 4");
    setDefault("corFilter",                    "expensive",  "Method to filter short reads from correction; 'quick' or 'expensive'; default 'expensive'");
    setDefault("corConsensus",                 "falcon",     "Which consensus algorithm to use; only 'falcon' is supported; default 'falcon'");

    #  Convert all the keys to lowercase, and remember the case-sensitive version

    foreach my $k (keys %synops) {
        (my $l = $k) =~ tr/A-Z/a-z/;

        $synnam{$l} = $k;                  #  Remember that option $l is stylized as $k.

        next  if (!exists($global{$k}));   #  If no option for this (it's a meta-option), skip.
        next  if ( exists($global{$l}));   #  If lowercase already exists, skip.

        $global{$l} = $global{$k};         #  Otherwise, set the lowercase option and
        delete $global{$k};                #  delete the uppercase version
    }

    #  If this is set, it breaks the consensus.sh and overlap.sh scripts.  Good grief!  Why
    #  are you running this in a task array!?

    if (exists($ENV{getGlobal("gridEngineTaskID")})) {
        undef $ENV{getGlobal("gridEngineTaskID")};
        print STDERR "ENV: ", getGlobal("gridEngineTaskID"), " needs to be unset, done.\n";
    }
}


#  Get the version information.  Needs to be last so that pathMap can be defined.

sub setVersion ($) {
    my $bin    = shift @_;
    my $version;

    open(F, "$bin/gatekeeperCreate --version 2>&1 |");
    while (<F>) {
        $version = $_;  chomp $version;
    }
    close(F);

    $global{'version'} = $version;
}


sub checkJava () {
    return  if ((getGlobal("corOverlapper") ne "mhap") &&
                (getGlobal("obtOverlapper") ne "mhap") &&
                (getGlobal("utgOverlapper") ne "mhap"));

    my $java       = getGlobal("java");
    my $versionStr = "unknown";
    my $version    = 0;

    #  Argh, we can't use runCommand() here, because we're included in Execution.pm.  Try to check
    #  it with -x.  Nope.  Fails if $java == "java".

    #if (! -x $java) {
    #    addCommandLineError("ERROR:  java executable '$java' not found or not executable\n");
    #}

    open(F, "$java -Xmx1g -showversion 2>&1 |");
    while (<F>) {
        #  First word is either "java" or "openjdk" or ...
        if (m/^.*\s+version\s+\"(\d+.\d+)(.*)\".*$/) {
            $versionStr = "$1$2";
            $version    =  $1;
        }
    }
    close(F);

    if ($version < 1.8) {
        addCommandLineError("ERROR:  mhap overlapper requires java version at least 1.8.0; you have $versionStr (from '$java').\n");
        addCommandLineError("ERROR:  '$java -showversion' reports:\n");

        open(F, "$java -Xmx1g -showversion 2>&1 |");
        while (<F>) {
            chomp;
            addCommandLineError("ERROR:    '$_'\n");
        }
        close(F);

    } else {
        print STDERR "-- Detected Java(TM) Runtime Environment '$versionStr' (from '$java').\n";
    }
}



sub checkGnuplot () {

    return  if (getGlobal("gnuPlotTested") == 1);

    my $gnuplot = getGlobal("gnuplot");
    my $format  = getGlobal("gnuplotImageFormat");
    my $version = undef;

    #  Check for existence of gnuplot.

    open(F, "$gnuplot -V |");
    while (<F>) {
        chomp;
        $version = $_;
        $version = $1  if ($version =~ m/^gnuplot\s+(.*)$/);
    }
    close(F);

    if (!defined($version)) {
        addCommandLineError("ERROR:  Failed to run gnuplot from '$gnuplot'.");
        addCommandLineError("ERROR:  Set option gnuplot=<path-to-gnuplot> or gnuplotTested=true to skip this test and not generate plots.\n");
        return;
    }

    #  Check for existence of a decent output format.  Need to redirect in /dev/null to make gnuplot
    #  not use it's builtin pager.

    if (!defined($format)) {
        my $havePNG = 0;
        my $haveSVG = 0;
        my $haveGIF = 0;

        open(F, "> /tmp/gnuplot-$$-test.gp");
        print F "set terminal\n";
        close(F);

        system("cd /tmp && $gnuplot < /dev/null /tmp/gnuplot-$$-test.gp > /tmp/gnuplot-$$-test.err 2>&1");

        open(F, "< /tmp/gnuplot-$$-test.err");
        while (<F>) {
            s/^\s+//;
            s/\s+$//;

            my @t = split '\s+', $_;

            $havePNG = 1  if ($t[0] eq 'png');
            $haveSVG = 1  if ($t[0] eq 'svg');
            $haveGIF = 1  if ($t[0] eq 'gif');
        }
        close(F);

        $format = "gif"   if ($haveGIF);
        $format = "svg"   if ($haveSVG);
        $format = "png"   if ($havePNG);

        setGlobal("gnuplotImageFormat", $format);

        unlink "/tmp/gnuplot-$$-test.gp";
        unlink "/tmp/gnuplot-$$-test.err";
    }

    if (!defined($format)) {
        addCommandLineError("ERROR:  Failed to detect a suitable output format for gnuplot.\n");
        addCommandLineError("ERROR:  Looked for png, svg and gif, found none of them.\n");
        addCommandLineError("Set option gnuplotImageFormat=<type>, or gnuplotTested=true to skip this test and not generate plots.\n");
        return;
    }

    #  Test if we can actually make images.

    open(F, "> /tmp/gnuplot-$$-test.gp");
    print F "set title 'gnuplot test'\n";
    print F "set xlabel 'X'\n";
    print F "set xlabel 'Y'\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output '/tmp/gnuplot-$$-test.1.$format'\n";
    print F "\n";
    print F "plot [-30:20] sin(x*20) * atan(x)\n\n";
    print F "\n";
    print F "set terminal $format size 256,256\n";
    print F "set output '/tmp/gnuplot-$$-test.2.$format'\n";
    print F "\n";
    print F "bogus line\n";
    close(F);

    #  Dang, we don't have runCommandSilently here, so have to do it the hard way.

    system("cd /tmp && $gnuplot < /dev/null /tmp/gnuplot-$$-test.gp > /tmp/gnuplot-$$-test.err 2>&1");

    if ((! -e "/tmp/gnuplot-$$-test.1.$format") ||
        (! -e "/tmp/gnuplot-$$-test.2.$format")) {
        addCommandLineError("ERROR:  gnuplot failed to generate images.\n");

        open(F, "< /tmp/gnuplot-$$-test.err");
        while (<F>) {
            chomp;
            addCommandLineError("ERROR:  gnuplot reports:  $_\n");
        }
        close(F);

        addCommandLineError("ERROR:  Set option gnuplotImageFormat=<type>, or gnuplotTested=true to skip this test and not generate plots.\n");
        return;
    }

    #  Yay, gnuplot works!

    print STDERR "-- Detected gnuplot version '$version' (from '$gnuplot') and image format '$format'.\n";
    #addCommandLineOption("gnuplotTested=1");

    unlink "/tmp/gnuplot-$$-test.gp";
    unlink "/tmp/gnuplot-$$-test.err";
    unlink "/tmp/gnuplot-$$-test.1.$format";
    unlink "/tmp/gnuplot-$$-test.2.$format";
}



sub checkParameters () {

    #
    #  Fiddle with filenames to make them absolute paths.
    #

    makeAbsolute("corOvlFrequentMers");
    makeAbsolute("obtOvlFrequentMers");
    makeAbsolute("utgOvlFrequentMers");

    #
    #  Adjust case on some of them
    #

    fixCase("corOverlapper");
    fixCase("obtOverlapper");
    fixCase("utgOverlapper");

    fixCase("corConsensus");
    fixCase("cnsConsensus");

    fixCase("corFilter");

    fixCase("unitigger");
    fixCase("stopAfter");

    #
    #  Well, crud.  'gridEngine' wants to be uppercase, not lowercase like fixCase() would do.
    #

    $global{"gridengine"} =~ tr/a-z/A-Z/;  #  NOTE: lowercase 'gridengine'

    #
    #  Check for inconsistent parameters
    #

    #  Genome size isn't properly decoded until later, but we want to fail quickly.  So, just test if
    #  a unitless number is supplied, and if that number is tiny.

    {
        my $gs = getGlobal("genomeSize");

        if (!defined($gs)) {
            addCommandLineError("ERROR:  Required parameter 'genomeSize' not set.\n");
        }

        if (($gs =~ m/^(\d+)$/) ||
            ($gs =~ m/^(\d+\.\d+)$/)) {
            if ($gs < 1000) {
                addCommandLineError("ERROR:  Implausibly small genome size $gs.  Check units!\n");
            }
        }
    }

    foreach my $var ("corOvlErrorRate", "obtOvlErrorRate", "utgOvlErrorRate", "corErrorRate", "obtErrorRate", "utgErrorRate", "cnsErrorRate") {
        if (!defined(getGlobal($var))) {
            addCommandLineError("ERROR:  Invalid '$var' specified; must be set\n");
        }
        elsif (getGlobal($var) !~ m/^[.-0123456789]/) {
            addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be numeric\n");
        }
        elsif ((getGlobal($var) < 0.0) || (getGlobal($var) > 1.0)) {
            addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be at least 0.0 and no more than 1.0\n");
        }
    }

    if (getGlobal("minReadLength") < getGlobal("minOverlapLength")) {
        my $mr = getGlobal("minReadLength");
        my $mo = getGlobal("minOverlapLength");

        addCommandLineError("ERROR:  minReadLength=$mr must be at least minOverlapLength=$mo.\n");
    }

    foreach my $var ("corOutCoverage") {
        if (!defined(getGlobal($var))) {
            addCommandLineError("ERROR:  Invalid 'corOutCoverage' specified; must be at least 1.0\n");
        }
        elsif (getGlobal($var) =~ m/all/i) {
            setGlobal($var, 9999);
        }
        elsif (getGlobal($var) !~ m/^[.-0123456789]/) {
            addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be numeric\n");
        }
        elsif (getGlobal($var) < 1.0) {
            addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be at least 1.0\n");
        }
    }

    foreach my $var ("corMaxEvidenceCoverageGlobal", "corMaxEvidenceCoverageLocal") {
        if (!defined(getGlobal($var))) {
            #  If undef, defaults to corOutCoverage in CorrectReads.pm
        }
        elsif (getGlobal($var) =~ m/^(\d*\.*\d*)(x*)$/) {
            if (($1 < 1.0) && ($2 ne "x")) {
                addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be at least 1.0\n");
            }
        }
        else {
            addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be numeric\n");
        }
    }

    foreach my $var ("utgGraphDeviation", "utgRepeatDeviation", "utgRepeatConfusedBP", "minReadLength", "minOverlapLength") {
        if (!defined(getGlobal($var))) {
            addCommandLineError("ERROR:  Invalid '$var' specified; must be set\n");
        }
        elsif (getGlobal($var) !~ m/^[.-0123456789]/) {
            addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be numeric\n");
        }
        elsif (getGlobal($var) < 0.0) {
            addCommandLineError("ERROR:  Invalid '$var' specified (" . getGlobal("$var") . "); must be at least 0.0\n");
        }
    }

    #
    #  Check for invalid usage
    #

    foreach my $tag ("cor", "obt", "utg") {
        if ((getGlobal("${tag}Overlapper") ne "mhap") &&
            (getGlobal("${tag}Overlapper") ne "ovl")  &&
            (getGlobal("${tag}Overlapper") ne "minimap")) {
            addCommandLineError("ERROR:  Invalid '${tag}Overlapper' specified (" . getGlobal("${tag}Overlapper") . "); must be 'mhap', 'ovl', or 'minimap'\n");
        }
    }

    foreach my $tag ("cor", "obt", "utg") {
        if (getGlobal("${tag}MhapSensitivity") eq "fast") {
            print STDERR "WARNING: deprecated ${tag}NhapSensitivity=fast replaced with ${tag}MhapSensitivity=low\n";
        }

        if ((getGlobal("${tag}MhapSensitivity") ne undef)    &&
            (getGlobal("${tag}MhapSensitivity") ne "low")    &&
            (getGlobal("${tag}MhapSensitivity") ne "normal") &&
            (getGlobal("${tag}MhapSensitivity") ne "high")) {
            addCommandLineError("ERROR:  Invalid '${tag}MhapSensitivity' specified (" . getGlobal("${tag}MhapSensitivity") . "); must be 'fast' or 'normal' or 'high'\n");
        }
    }

    if ((getGlobal("ovsMethod") ne "sequential") &&
        (getGlobal("ovsMethod") ne "parallel")) {
        addCommandLineError("ERROR:  Invalid 'ovsMethod' specified (" . getGlobal("ovsMethod") . "); must be 'sequential' or 'parallel'\n");
    }
    if ((getGlobal("useGrid")   eq "0") &&
        (getGlobal("ovsMethod") eq "parallel")) {
        addCommandLineError("ERROR:  ovsMethod=parallel requires useGrid=true or useGrid=remote.  Set ovsMethod=sequential if no grid is available\n");
    }

    if ((getGlobal("unitigger") ne "bogart")) {
        addCommandLineError("ERROR:  Invalid 'unitigger' specified (" . getGlobal("unitigger") . "); must be 'unitigger' or 'bogart'\n");
    }

    if ((getGlobal("corConsensus") ne "falcon")) {
        addCommandLineError("ERROR:  Invalid 'corConsensus' specified (" . getGlobal("corConsensus") . "); must be 'utgcns' or 'falcon' or 'falconpipe'\n");
    }

    if ((getGlobal("cnsConsensus") ne "quick") &&
        (getGlobal("cnsConsensus") ne "pbdagcon") &&
        (getGlobal("cnsConsensus") ne "utgcns")) {
        addCommandLineError("ERROR:  Invalid 'cnsConsensus' specified (" . getGlobal("cnsConsensus") . "); must be 'quick', 'pbdagcon', or 'utgcns'\n");
    }


    if ((!defined("lowCoverageAllowed") &&  defined("lowCoverageDepth")) ||
        ( defined("lowCoverageAllowed") && !defined("lowCoverageDepth"))) {
        addCommandLineError("ERROR:  Invalid 'lowCoverageAllowed' and 'lowCoverageDepth' specified; both must be set\n");
    }

    if ((getGlobal("saveOverlaps") ne "0") &&
        (getGlobal("saveOverlaps") ne "stores") &&
        (getGlobal("saveOverlaps") ne "1")) {
        addCommandLineError("ERROR:  Invalid 'saveOverlaps' specified (" . getGlobal("saveOverlaps") . "); must be 'false', 'stores', or 'true'\n");
    }

    if ((getGlobal("corFilter") ne "quick") &&
        (getGlobal("corFilter") ne "expensive") &&
        (getGlobal("corFilter") ne "none")) {
        addCommandLineError("ERROR:  Invalid 'corFilter' specified (" . getGlobal("corFilter") . "); must be 'none' or 'quick' or 'expensive'\n");
    }


    if ((getGlobal("useGrid") ne "0") &&
        (getGlobal("useGrid") ne "1") &&
        (getGlobal("useGrid") ne "remote")) {
        addCommandLineError("ERROR:  Invalid 'useGrid' specified (" . getGlobal("useGrid") . "); must be 'true', 'false' or 'remote'\n");
    }

    if (defined(getGlobal("stopAfter"))) {
        my $ok = 0;
        my $st = getGlobal("stopAfter");
        $st =~ tr/A-Z/a-z/;

        my $failureString = "ERROR:  Invalid stopAfter specified (" . getGlobal("stopAfter") . "); must be one of:\n";

        my @stopAfter = ("gatekeeper",
                         "meryl",
                         "overlapConfigure",
                         "overlap",
                         "overlapStoreConfigure",
                         "overlapStore",
                         "readCorrection",
                         "readTrimming",
                         "unitig",
                         "consensusConfigure",
                         "consensusCheck",
                         "consensusLoad",
                         "consensusAnalyze");

        foreach my $sa (@stopAfter) {
            $failureString .= "ERROR:      '$sa'\n";
            $sa =~ tr/A-Z/a-z/;
            if ($st eq $sa) {
                $ok++;
                setGlobal('stopAfter', $st);
            }
        }

        addCommandLineError($failureString)   if ($ok == 0);
    }

    {
        my @v = split '\s+', getGlobal("contigFilter");

        if (scalar(@v) != 5) {
            addCommandLineError("contigFilter must have five values: minReads minLength singleReadSpan lowCovFraction lowCovDepth\n");
        }

        addCommandLineError("contigFilter 'minReads' must be a positive integer, currently $v[0]\n")          if (($v[0] < 0) || ($v[0] !~ m/^[0-9]+$/));
        addCommandLineError("contigFilter 'minLength' must be a positive integer, currently $v[1]\n")         if (($v[1] < 0) || ($v[1] !~ m/^[0-9]+$/));
        addCommandLineError("contigFilter 'singleReadSpan' must be between 0.0 and 1.0, currently $v[2]\n")   if (($v[2] < 0) || (1 < $v[2]) || ($v[2] !~ m/^[0-9]*\.{0,1}[0-9]*$/));
        addCommandLineError("contigFilter 'lowCovFraction' must be between 0.0 and 1.0, currently $v[3]\n")   if (($v[3] < 0) || (1 < $v[3]) || ($v[3] !~ m/^[0-9]*\.{0,1}[0-9]*$/));
        addCommandLineError("contigFilter 'lowCovDepth' must be a positive integer, currently $v[4]\n")       if (($v[4] < 0) || ($v[4] !~ m/^[0-9]+$/));
    }
}


1;
