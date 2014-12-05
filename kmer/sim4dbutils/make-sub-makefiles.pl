#!/usr/bin/perl

use strict;

sub makemake ($@) {
    my $f = shift @_;

    open(F, "> $f.mk");
    print F "TARGET  := $f\n";
    print F "SOURCES := $f.C\n"    if (scalar(@_) == 0);
    print F "SOURCES := $f.C @_\n" if (scalar(@_) > 0);
    print F "\n";
    print F "SRC_INCDIRS := ../libutil ../libbio ../libseq ../libsim4\n";
    print F "\n";
    print F "TGT_LDFLAGS := -L\${TARGET_DIR}\n";
    print F "TGT_LDLIBS  := -lseq -lsim4 -lbio -lutil\n";
    print F "TGT_PREREQS := libseq.a libbio.a libutil.a libsim4.a\n";
    close(F);

    print STDOUT "$f.mk\n";
}

makemake("comparePolishes", "s4p_overlap.C");
makemake("removeDuplicate", "s4p_overlap.C");
makemake("removeRedundant", "s4p_overlap.C");

makemake("cleanPolishes");
makemake("fixPolishesIID");
makemake("convertToAtac");
makemake("convertToExtent");
makemake("convertPolishes");
makemake("detectChimera");
makemake("depthOfPolishes");
makemake("filterPolishes");
makemake("headPolishes");
makemake("mappedCoverage");
makemake("mergePolishes");
makemake("parseSNP");
makemake("pickBestPolish");
makemake("pickBestPair");
makemake("pickUniquePolish");
makemake("plotCoverageVsIdentity");
makemake("sortPolishes");
makemake("summarizePolishes");
makemake("uniqPolishes");
makemake("vennPolishes");
makemake("realignPolishes");
makemake("reportAlignmentDifferences");

