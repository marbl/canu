#!/usr/bin/perl

use strict;

$| = 1;

my $leaff = "./leaff";

#
#  A leaff test, using a large file.
#

#  Make a big script to write a large random file, using unique defline names
#
my @names = ("zer",
             "one",
             "two",
             "thr",
             "fou",
             "fiv",
             "six",
             "sev",
             "eig",
             "nin");

if (! -e "leaff-test-large.scr") {
    print STDERR "Writing a script for leaff\n";

    open(F, "> leaff-test-large.scr");
    for (my $i="000000"; $i ne "1000000"; $i++) {
        my $id = $i;

        $id =~ s/0/zer/g;
        $id =~ s/1/one/g;
        $id =~ s/2/two/g;
        $id =~ s/3/thr/g;
        $id =~ s/4/fou/g;
        $id =~ s/5/fiv/g;
        $id =~ s/6/six/g;
        $id =~ s/7/sev/g;
        $id =~ s/8/eig/g;
        $id =~ s/9/nin/g;

        print F "-h $id -G 1 500 300\n";
    }
    close(F);
    print STDERR "\n";
}

 
if (! -e "leaff-test-large.fasta") {
    print STDERR "Running leaff to build fasta-file\n";
    system("${leaff} -A leaff-test-large.scr > leaff-test-large.fasta");
}

#if (! -e "leaff-test-large.fastaidx") {
#    print STDERR "Running leaff to build iid index\n";
#    system("${leaff} -F leaff-test-large.fasta");
#}

if (! -e "leaff-test-large.fastaidx") {
    print STDERR "Running leaff to build defline index\n";
    system("${leaff} -Fd leaff-test-large.fasta");
}




#  Pull out some sequences
#
print STDERR "Finding zerzerzeroneoneone\n";
system("${leaff} -Fd leaff-test-large.fasta -s zerzerzeroneoneone");

print STDERR "Finding ninninninzerzerzer to ninninninzerzertwo\n";
system("${leaff} -Fd leaff-test-large.fasta -S ninninninzerzerzer ninninninzerzertwo");

print STDERR "Finding ninninninninninnin (the last sequence)\n";
system("${leaff} -Fd leaff-test-large.fasta -s ninninninninninnin");

print STDERR "Finding zerzerzerzerzerzer (the first sequence)\n";
system("${leaff} -Fd leaff-test-large.fasta -s zerzerzerzerzerzer");

print STDERR "Finding last, then first, using iid's\n";
system("${leaff} -Fd leaff-test-large.fasta -Ii -s 999999 -s 0");


