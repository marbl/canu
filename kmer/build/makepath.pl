#!/usr/bin/env perl

use Cwd 'abs_path';

if ($#ARGV<0) {
    printf STDERR "usage: makepath.pl path [subpath]\n";
    exit(1);
}

@path = split('/',abs_path($ARGV[0]));
if (!@path) {
    print STDERR "makepath.pl: The path '",$ARGV[0],"' does not exist.\n";
    print STDERR "makepath.pl: perhaps you are building from an incomplete CVS?\n";
    print "\n";
    exit(0);
}

@subpath = ();
if ($#ARGV) {
    @subpath = split('/',abs_path($ARGV[1])) ;
    $subpath[0] = '';
}

while(@subpath && @path && $path[0] eq $subpath[0]) {
    shift(@path); 
    shift(@subpath);
}
for (@subpath) {$_ = '..'};

print join("/", @subpath, @path), "/\n";
