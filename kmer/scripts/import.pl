#!/bin/perl

#  Script that will import a genome into the GENOMES directory.
#
#  usage:
#    import.pl NICKNAME /path/to/genome.fasta
#

use strict;

my $leaff = "/usr/local/ir/bin/leaff";

if (scalar(@ARGV) != 2) {
    print STDERR "usage: $0 NICKNAME genome.fasta\n";
    exit(1);
}

my $nickname = shift @ARGV;
my $genome   = shift @ARGV;

if (-e "/IR/devel/genomics/GENOMES/$nickname") {
    print STDERR "error: /IR/devel/genomics/GENOMES/$nickname exists.\n";
    exit(1);
}

system("mkdir /IR/devel/genomics/GENOMES/$nickname");

my $d = `date`;    chomp $d;
my $u = `whoami`;  chomp $u;

open(F, "> /IR/devel/genomics/GENOMES/$nickname/$nickname.readme");
print F "species     = \n";
print F "type        = \n";
print F "source      = \n";
print F "version     = \n";
print F "releasedate = \n";
print F "mappingdate = \n";
print F "notes:\n";
print F "Original source file is $genome.\n";
print F "Imported on $d by $u.\n";
close(F);


#  Squeeze the file, make some extra files
#
#print "$leaff -f $genome -W > /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta\n";
system("$leaff -f $genome -W > /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta");

#print "$leaff -Fdc /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta\n";
system("$leaff -Fdc /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta");

#print "grep '>' /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta > /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta.deflines\n";
system("grep '>' /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta > /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta.deflines");

#print "/work/assembly/walenzbp/releases/dumpBlocks /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta > /IR/devel/genomics/GENOMES/$nickname/$nickname.breakpoints\n";
system("/work/assembly/walenzbp/releases/dumpBlocks /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta > /IR/devel/genomics/GENOMES/$nickname/$nickname.breakpoints");



#  Write out the atai index
#
#print "$leaff -Fd /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta -i $nickname\n";

if (! -e "/IR/devel/genomics/GENOMES/$nickname/$nickname.atai") {
    open(O, "> /IR/devel/genomics/GENOMES/$nickname/$nickname.atai");
    open(F, "$leaff -Fd /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta -i $nickname |");

    #  Format line and S record
    #
    my $junk;
    $junk = <F>;  print O $junk;
    $junk = <F>;  print O $junk;

    while (<F>) {
        chomp;

        if (m/^G\s+\w+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(.*)$/) {
            my $a    = "G ch $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 $11";
            my $name = $12;
            my $uid  = $13;
            my $def  = $14;

            #print "'$12' '$13' '$14'\n";

            if (length($def) > 0) {
                if      ($def =~ m/chromosome=(Chr\S+)/) {
                    $name = $1;
                } elsif ($def =~ m/chromosome=(\S+)/) {
                    $name = "Chr$1";
                } elsif ($def =~ m/GA_(x\S+):/) {
                    $name = $1;
                } elsif ($def =~ m/^>chr(.*)$/) {
                    $name = "Chr$1";
                } elsif ($def =~ m/(NT.\d+.\d+)/) {
                    $name = $1;
                } elsif ($def =~ m/ga_uid=(\d+)/) {
                    $name = $1;
                }


                if      ($def =~ m/ga_uid=(\d+)/) {
                    $uid = $1;
                } elsif ($def =~ m/^>(\d+)/) {
                    $uid = $1;
                } elsif ($def =~ m/^>CRA\|(\S+)/) {
                    $uid = $1;
                } elsif ($def =~ m/gi\|(\d+)/) {
                    $uid = $1;
                } elsif ($def =~ m/(Rnor.*)\|/) {
                    $uid = $1;
                }
            }

            #print "$name -- $uid\n";

            print O "$a $name $uid $def\n";
        } else {
            print STDERR "WARNING:  No match '$_'.\n";
            print O "$_\n";
        }
    }
    close(F);
    close(O);
}


print STDERR "Edit /IR/devel/genomics/GENOMES/$nickname/$nickname.readme\n";
