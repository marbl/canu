#!/usr/bin/env perl

use Cwd 'abs_path';

if ($#ARGV<0) {
    printf STDERR "usage: makepath.pl path [subpath]\n";
    exit(1);
}

@path = split('/',abs_path($ARGV[0]));
@subpath = ();
if ($#ARGV) {
    @subpath = split('/',abs_path($ARGV[1])) ;
    $subpath[0] = '';
}

while(@subpath && @path && $path[0] eq $subpath[0]) {
    shift(@path); 
    shift(@subpath);
}

print "../" x @subpath, join("/", @path), "/\n";
