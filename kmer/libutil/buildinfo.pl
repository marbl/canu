#!/usr/bin/perl

use strict;

my $functionName = shift @ARGV;
my $compileOpts  = shift @ARGV;

print STDERR "Creating build information with name '$functionName'.\n";

open(Z, "> buildinfo-$functionName.h");
print Z "#include <stdio.h>\n";
print Z "#ifdef __cplusplus\n";
print Z "extern \"C\" {\n";
print Z "#endif\n";
print Z "void\n";
print Z "buildinfo_$functionName(FILE *F);\n";
print Z "#ifdef __cplusplus\n";
print Z "}\n";
print Z "#endif\n";
close(Z);

open(Z, "> buildinfo-$functionName.c");

print Z "#include <stdio.h>\n";
print Z "\n";
print Z "void\n";
print Z "buildinfo_$functionName(FILE *F) {\n";
print Z "  fprintf(F, \"===================================================================\\n\");\n";
print Z "  fprintf(F, \"$functionName\\n\");\n";
print Z "  fprintf(F, \"  built on %s %s.\\n\", __DATE__, __TIME__);\n";
print Z "  fprintf(F, \"  built with $compileOpts\\n\");\n";
print Z "  fprintf(F, \"\\n\");\n";

my @status;

foreach my $f (@ARGV) {
    if ($f !~ m/buildinfo_/) {
        open(F, "cvs status -v $f |");
        while (<F>) {
            chomp;
            print Z "  fprintf(F, \"$_\\n\");\n" if (! m/^$/) && (! m/^=+$/);
        }
    }
}

print Z "}\n";

close(Z);
