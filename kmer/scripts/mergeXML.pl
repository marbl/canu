#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

#
#  Merges directories of XML into one.  The layout of the directories
#  is irrelevant, but the names of the files is relevant.
#
#  It will find all the files in the directories, and merge
#  those that have the same scaffold id:
#
#  /some/path/to/the/directory/namespace:uid.gbf
#
#  The scaffold id is found with the regex /(\d+).gbf$/
#


if (scalar @ARGV < 1) {
    print STDOUT "usage: $0 outdir indir indir indir ...\n";
    exit;
}

my $outdir = shift @ARGV;
system("mkdir $outdir") if (! -d $outdir);

my %inputfiles;

#
#  Find all the input files, and store them sorted by scaffold.
#
foreach my $i (@ARGV) {
    open(F, "find $i -type f -name '*.gbf' -print |");
    while (<F>) {
        chomp;
        if (m/(\d+).gbf/) {
            if (defined($inputfiles{$1})) {
                $inputfiles{$1} .= "\0$_";
            } else {
                $inputfiles{$1} = "$_";
            }
        } else {
            print STDERR "Invalid filename: '$_'\n";
        }
    }
    close(F);
}


#
#  For each scaffold we found, write new output
#
foreach my $k (keys %inputfiles) {
    my @files = split '\0', $inputfiles{$k};

    #print STDERR "Merging $k with ", scalar @files, " files.\n";

    #  Check that all the inputs are the same assembly
    #
    my %check;
    my $headerLine;
    foreach my $f (@files) {
        open(F, "< $f");
        $headerLine = <F>;
        chomp $headerLine;
        $check{$headerLine} = 1;
    }

    if (scalar(keys %check) > 1) {
        print STDERR "ERROR:  Inconsistency in headers for $k:\n";
        foreach my $x (keys %check) {
            print STDERR "        $x\n";
        }
        exit;
    }


    $k =~ m/(\d\d)$/;
    my $dir = $1;

    system("mkdir $outdir/$dir") if (! -d "$outdir/$dir");

    open(Z, "> $outdir/$dir/$k.gbf");
    print Z "$headerLine\n";

    foreach my $f (@files) {
        open(F, "< $f");
        $_ = <F>;
        while (<F>) {
            if ($_ ne "</game>\n") {
                print Z $_;
            }
        }
        close(F);
    }

    print Z "</game>\n";
    close(Z);
}
