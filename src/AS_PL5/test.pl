#! /usr/bin/perl

use strict;
use warnings;
use AS::MSG;
use Devel::Peek;

my $parser = AS::MSG::Parser->new();
my $i = 0;
my @filters = ("MDI", "CLK");
$parser->parser_initialize(*STDIN{IO});
$parser->{"filter_out"} = \@filters;

@filters = @{$parser->{"filter_in"}};
print("in: @filters\n");
@filters = @{$parser->{"filter_out"}};
print("out: @filters\n");

while($parser->parser_next()) {
    $i++;
}

print($i, "\n");
