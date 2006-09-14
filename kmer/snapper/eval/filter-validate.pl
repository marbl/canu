#!/bin/perl

#
#  Read in the output of a snapper2 validation run, print the highest
#  specificity for each distinct sensitivity.
#

my $numToShow = 10;
if (! -e $ARGV[0]) {
    $numToShow = shift @ARGV;
}

my $hdrshown = 0;

foreach my $file (@ARGV) {
    my %spec;
    my %line;

    open(F, "< $file");
    my $hdr = <F>;
    chomp $hdr;

    while (<F>) {
        chomp;

        my @vals = split '\s+', $_;

        #  3 -> sensitivity
        #  4 -> specificity

        if ($spec{$vals[3]} < $vals[4]) {
            $spec{$vals[3]} = $vals[4];
            $line{$vals[3]} = $_;
        }
    }

    close(F);

#        print "\n$file\n                 $hdr\n";
#        my @sortedK = sort { $b <=> $a } keys %spec;
#        $#sortedK = $numToShow - 1;
#        foreach my $k (@sortedK) {
#            print "$k $spec{$k} -- $line{$k}\n";
#        }

    if ($hdrshown == 0) {
        print "                                                         $hdr\n";
        $hdrshown = 1;
    }

    $file = substr("$file                              ", 0, 40);
    my @sortedK = sort { $b <=> $a } keys %spec;
    $#sortedK = $numToShow - 1;
    foreach my $k (@sortedK) {
        printf "$file$k $spec{$k} -- $line{$k}\n";
    }

    undef @sortedK;
    undef %spec;
}
