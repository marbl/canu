#!/bin/perl

#  Import a genome into the IR GENOMES directory.
#
#  usage:
#    import.pl NICKNAME /path/to/genome.fasta
#
#  Revisions:
#    20040419 - bpw - check the results of system() calls.
#

use strict;

my $leaff = "/usr/local/ir/bin/leaff";

if (scalar(@ARGV) != 2) {
    print STDERR "usage: $0 NICKNAME genome.fasta\n";
    exit(1);
}
my $nickname = shift @ARGV;
my $genome   = shift @ARGV;
my $d        = `date`;    chomp $d;
my $u        = `whoami`;  chomp $u;

if (-e "/IR/devel/genomics/GENOMES/$nickname") {
    print STDERR "error: /IR/devel/genomics/GENOMES/$nickname exists.\n";
    exit(1);
}

if (! -x $leaff) {
    print STDERR "Can't find or execute $leaff\n";
    exit(1);
}

system("mkdir /IR/devel/genomics/GENOMES/$nickname") and die "Failed to create /IR/devel/genomics/GENOMES/$nickname";
system("chmod 775 /IR/devel/genomics/GENOMES/$nickname") and die "Failed to chmod\n";

open(F, "> /IR/devel/genomics/GENOMES/$nickname/$nickname.readme") or die "Failed to open readme\n";
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
system("$leaff -f $genome -W > /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta") and die "Failed to squeeze the genome.\n";
system("$leaff -Fd /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta") and die "Failed to index the genome.\n";
system("grep '>' /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta > /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta.deflines") and die "Failed to find deflines.\n";
system("/usr/local/ir/bin//dumpBlocks /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta > /IR/devel/genomics/GENOMES/$nickname/$nickname.breakpoints") and die "Failed to get sequence blocks.\n";


#  Write out the atai index
#
#print "$leaff -Fd /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta -i $nickname\n";

if (! -e "/IR/devel/genomics/GENOMES/$nickname/$nickname.atai") {
    open(O, "> /IR/devel/genomics/GENOMES/$nickname/$nickname.atai") or die "Failed to open atai file for output.\n";
    open(F, "$leaff -Fd /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta -i $nickname |") or die "Failed to dump the genome index.\n";

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

open(F, ">> assemblies.atai") and die "Failed to append the info line to assemblies.atai\n";
print F "S $nickname /IR/devel/genomics/GENOMES/$nickname/$nickname.fasta\n";
close(F);

print STDERR "Please fill out the information in /IR/devel/genomics/GENOMES/$nickname/$nickname.readme\n";
