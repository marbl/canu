#!/usr/bin/perl

# Handle the -o arguement for msvc...

@NEWARGS = ();

$cwd = `pwd`;
($cwdhead = $cwd) =~ s/[\/\\][^\/\\]*$//;

for ($i=0;$i<@ARGV;$i++) {
    if ($ARGV[$i] eq "-o") {
	$dest = $ARGV[$i+1];
	die if ($dest eq "");
	($path = $dest) =~ s/[\/\\][^\/\\]*$//;
	($base = $dest) =~ s/^.*[\/\\]([^\/\\]*)$/\1/;
	($basenoext = $base) =~ s/[.][^.]*$//;
	$i++;
	push @NEWARGS, "-Fo$dest";
    } else {
	push @NEWARGS, "$ARGV[$i]";
    }
}

$cmd = "cl";
for ($i=0;$i<@NEWARGS;$i++) {
    $cmd .= " $NEWARGS[$i]";
}
print "$cmd","\n";
exec "$cmd";
