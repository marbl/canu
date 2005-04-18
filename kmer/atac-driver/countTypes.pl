#!/usr/bin/perl

my $leng;
my $last;
my $this;
my $type;

my $bpnew;
my $bpext;
my $bpdiffext;
my $bpunk;

my $count = 0;

open(F, "< overlap.Aannotation");
while (<F>) {
    chomp;
    if (m/(.)\s+\d+\s*:\s*\d+\s*-\s*\d+\s*\[\s*(\d+)\]\s/) {
        $last = $this;
        $this = $1;
        $leng = $2;
    } else {
        print STDERR "Failed on: '$_'\n";
    }

    if      (($last eq "U") && ($this eq "1")) {
        $type = "new";
    } elsif ($this eq "U") {
        undef $type;
    } elsif (($last eq "Y") && ($this eq "1")) {
        $type = "ext";
    } elsif ($this eq "1") {
        $type = "diffext";
    } elsif (($last eq "Y") && ($this eq "Y")) {
    } elsif (($last eq "N") && ($this eq "N")) {
    } else {
        $type = "unknown";
    }

    if ($type eq "new") {
        $bpnew += $leng;
    } elsif ($type eq "ext") {
        $bpext += $leng;
    } elsif ($type eq "diffext") {
        $bpdiffext += $leng;
    } elsif ($type eq "unknown") {
        $bpunk += $leng;
    }

    $count++;
    if (($count % 100000) == 0) {
        print "$count\n";
        print "NEW bp:              $bpnew\n";
        print "EXTENSION            $bpext\n";
        print "DIFFERENT EXTENSION  $bpdiffext\n";
        print "unkown               $bpunk\n";
    }
}
close(F);

print "NEW bp:              $bpnew\n";
print "EXTENSION            $bpext\n";
print "DIFFERENT EXTENSION  $bpdiffext\n";
print "unkown               $bpunk\n";
