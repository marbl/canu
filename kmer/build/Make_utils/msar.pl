#!/usr/bin/perl

# Handle the output argument for msvc and lib building...

@NEWARGS = ();

$cwd = `pwd`;
($cwdhead = $cwd) =~ s/[\/\\][^\/\\]*$//;

$seenlib = 0;
for ($i=0;$i<@ARGV;$i++) {
    if (!$seenlib && $ARGV[$i] !~ /^-/ && $ARGV[$i] =~ /\.lib$/ ) {
	$dest = $ARGV[$i];
	die if ($dest eq "");
	($path = $dest) =~ s/[\/\\][^\/\\]*$//;
	($base = $dest) =~ s/^.*[\/\\]([^\/\\]*)$/\1/;
	($basenoext = $base) =~ s/[.][^.]*$//;
	push @NEWARGS, "-out:$dest";
	$seenlib = 1;
    } else {
	push @NEWARGS, "$ARGV[$i]";
    }
}

$cmd = "lib";
for ($i=0;$i<@NEWARGS;$i++) {
    $cmd .= " $NEWARGS[$i]";
}
print "$cmd","\n";
exec "$cmd";
