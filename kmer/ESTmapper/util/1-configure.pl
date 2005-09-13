use strict;


#  Takes three args:
#    path to our working directory
#    path to the genome build directory
#    path to the ESTs
#
sub configure {
    my $path      = shift @_;
    my $genomic   = shift @_;
    my $cdna      = shift @_;

    print STDERR "ESTmapper: Performing a configure.\n";

    ($path eq "")    and die "ERROR: ESTmapper/configure-- no directory given.\n";
    ($genomic eq "") and die "ERROR: ESTmapper/configure-- no genomic sequences given.\n";
    ($cdna eq "")    and die "ERROR: ESTmapper/configure-- no cDNA sequences given.\n";
    (! -f $cdna)     and die "ERROR: ESTmapper/configure-- can't find the cdna sequence '$cdna'\n";

    #  Make some organization
    #
    mkdir "$path"          if (! -d "$path");
    mkdir "$path/0-input"  if (! -d "$path/0-input");
    mkdir "$path/1-search" if (! -d "$path/1-search");
    mkdir "$path/2-filter" if (! -d "$path/2-filter");
    mkdir "$path/3-polish" if (! -d "$path/3-polish");

    #  XXX:  We should check that the genome dir is valid and complete.
    #
    system("ln -s $genomic $path/0-input/genome");

    #  Check the input files exist, create symlinks to them, and find/build index files
    #
    symlink "${cdna}",       "$path/0-input/cDNA.fasta"       if ((! -f "$path/0-input/cDNA.fasta"));
    symlink "${cdna}idx",    "$path/0-input/cDNA.fastaidx"    if ((! -f "$path/0-input/cDNA.fastaidx") && (-f "${cdna}idx"));

    if (! -f "$path/0-input/cDNA.fastaidx") {
        print STDERR "ESTmapper/configure-- Generating the index for '$path/0-input/cDNA.fasta'\n";
        runCommand("$leaff -F $path/0-input/cDNA.fasta") and die "Failed.\n";
    }
}

1;
