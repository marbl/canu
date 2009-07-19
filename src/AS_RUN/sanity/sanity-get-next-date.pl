#!/usr/bin/perl

use Time::Local;

my @v;

my $date = shift @ARGV;
my $offt = shift @ARGV;
my $type = shift @ARGV;
my $time = time();

if ($date =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
    $time = timelocal(0, $5, $4, $3, $2 - 1, $1);
}

$time += $offt;

my @v = localtime($time);

$v[5] += 1900;
$v[4]++;

$v[5] = substr("0000$v[5]", -4);
$v[4] = substr("0000$v[4]", -2);
$v[3] = substr("0000$v[3]", -2);
$v[2] = substr("0000$v[2]", -2);
$v[1] = substr("0000$v[1]", -2);

#$v[2] = "00";
#$v[1] = "01";

if      ($type eq "next") {
    $thisdate = "$v[5]-$v[4]-$v[3]-$v[2]$v[1]";
} elsif ($type eq "hold") {
    $thisdate = "$v[5]$v[4]$v[3]$v[2]$v[1].00";
} else {
    $thisdate = "$v[4]$v[3]-$v[2]$v[1]";
}

print "$thisdate\n";


