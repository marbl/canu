#!/usr/bin/env perl

use strict;

if (scalar(@ARGV) < 3) {
    print STDERR "usage: $0 functionName compilerOptions sourceFile sourceFile ...\n";
    exit(1);
}

my $functionName = shift @ARGV;
my $compileOpts  = shift @ARGV;

#print STDERR "Creating build information with name '$functionName'.\n";

open(Z, "> buildinfo-$functionName.h") or die "Failed to open buildinfo-$functionName.h\n";
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

open(Z, "> buildinfo-$functionName.c") or die "Failed to open buildinfo-$functionName.c\n";;
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

#  Check that there is a CVS root; if not found, assume we are
#  building a copied tree.
#
open(F, "< CVS/Root");
my $root = <F>;
chomp $root;
close(F);

if (-d $root) {
    foreach my $f (@ARGV) {
        if ($f !~ m/buildinfo_/) {
            open(F, "cvs status -v $f |");
            while (<F>) {
                chomp;
                print Z "  fprintf(F, \"$_\\n\");\n" if (! m/^$/) && (! m/^=+$/);
            }
        }
    }
} else {
    print Z "  fprintf(F, \"  No source repository found during build; status not available\\n\");\n";
}

print Z "}\n";

close(Z);
