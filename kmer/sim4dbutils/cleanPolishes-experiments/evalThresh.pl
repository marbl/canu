#!/bin/perl

$threshold = 10000;

while ($threshold < 500000) {
    print "THRESHOLD = $threshold\n";
    system("./splitMatches -qquiet -threshold $threshold < /part3/polishes-good");
    $threshold += 10000
}
