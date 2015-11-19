#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use MIME::QuotedPrint;
use MIME::Base64;
use HTML::Entities;
use strict;

if (scalar(@ARGV) != 5) {
    die "wrong args.\n";
}

my $wrkdir     = shift @ARGV;
my $thisdate   = shift @ARGV;
my $label      = shift @ARGV;
my $addresses  = shift @ARGV;
my $assemblies = shift @ARGV;
my @assemblies = split ',', $assemblies;
my $timenow    = localtime();

print "Subject: CAtest $thisdate $label\n";
print "To:      $addresses\n";
print "Content-Type: multipart/mixed; boundary=Bri_Says_This_Is_The_Boundary\n";  # alternative if no atachments
print "\n";


########################################
#
#  Dump the HTML part
#
print "--Bri_Says_This_Is_The_Boundary\n";
print "Content-Type: text/html; charset=\"iso-8859-1\"\n";
print "Content-Transfer-Encoding: quoted-printable\n";
print "\n";

print "<HTML>\n";
print "<P>\n";
print "<font size=4>" . encode_entities("Results for $thisdate (Finished at $timenow).") . "</font>\n";
print "</P>\n";
print "\n";

print "<P>\n";

foreach my $asm (@assemblies) {
    open(R, "< $wrkdir/$thisdate/$asm/sanity-result.out");
    $_ = <R>;
    close(R);

    chomp;
    #  Green!
    s!SUCCESS!<span style=\"color: rgb(0, 255, 0);\">SUCCESS</span>!g;
    s!same!<span style=\"color: rgb(0, 255, 0);\">same</span>!g;

    #  Red
    s!FAILURE!<span style=\"color: rgb(255, 0, 0);\">FAILURE</span>!g;
    s!differs!<span style=\"color: rgb(255, 0, 0);\">differs</span>!g;

    print "$_<BR>\n";
}

print "</P>\n";

print "<HR>\n";
if ((! -e "$wrkdir/$thisdate/wgs/kmer.updates") || (-z "$wrkdir/$thisdate/wgs/kmer.updates")) {
    print "<P>\n";
    print "<font size=4>No changes to kmer.</font>\n";
    print "</P>\n";
} else {
    print "<P>\n";
    print "<font size=4>Changes for kmer.</font>\n";
    print "</P>\n";
    print "<P>\n";
    open(F, "< $wrkdir/$thisdate/wgs/kmer.updates");
    while (<F>) {
        chomp;
        print encode_entities($_) . "<BR>\n";
    }
    close(F);
    print "</P>\n";
}

print "<HR>\n";
if ((! -e "$wrkdir/$thisdate/wgs/src.updates") || (-z "$wrkdir/$thisdate/wgs/src.updates")) {
    print "<P>\n";
    print "<font size=+2>No chages to wgs-assembler.</font>\n";
    print "</P>\n";
} else {
    print "<P>\n";
    print "<font size=+2>Changes for wgs-assembler.</font>\n";
    print "</P>\n";
    print "<P>\n";
    open(F, "< $wrkdir/$thisdate/wgs/src.updates");
    while (<F>) {
        chomp;
        print encode_entities($_) . "<BR>\n";
    }
    close(F);
    print "</P>\n";
}

#print "<HR>\n";
#print "<P>\n";
#print "<font size=+2>Build Results for kmer.</font>\n";
#print "</P>\n";
#open(F, "< $wrkdir/$thisdate/wgs/kmer/make.err");
#while (<F>) {
#    chomp;
#    print encode_entities($_) . "<BR>\n";
#}
#close(F);

#print "<HR>\n";
#print "<P>\n";
#print "<font size=+2>Build Results for wgs-assembler.</font>\n";
#print "</P>\n";
#open(F, "< $wrkdir/$thisdate/wgs/src/make.err");
#while (<F>) {
#    chomp;
#    print encode_entities($_) . "<BR>\n";
#}
#close(F);

print "</HTML>\n";
print "\n";



########################################
#
#  Dump the ASCII part -- not dumped; outlook shows this as another attachment
#
#print "--Bri_Says_This_Is_The_Boundary\n";
#print "Content-Type: text/plain; charset=\"iso-8859-1\"\n";
#print "Content-Transfer-Encoding: quoted-printable\n";
#print "\n";



########################################
#
#  Dump the checkout and build results
#
attachFile("$wrkdir/$thisdate/wgs/kmer.updates",  "kmerChanges.txt");
attachFile("$wrkdir/$thisdate/wgs/src.updates",   "wgsChanges.txt");
attachFile("$wrkdir/$thisdate/wgs/kmer/make.err", "kmerBuildError.txt");
attachFile("$wrkdir/$thisdate/wgs/src/make.err",  "wgsBuildError.txt");


########################################
#
#  Dump the QC attachments
#
foreach my $asm (@assemblies) {
    my $errfile = "$wrkdir/$thisdate/$asm/sanity-error.out";
    my $qcfile  = "$wrkdir/$thisdate/$asm/sanity-qc.out";

    if ((-e $qcfile) && (! -z $qcfile)) {
        my @n       = split "/", $qcfile;
        my $l       = scalar(@n);
        my $name    = "$n[$l-3]-$n[$l-2]-qc.txt";

        attachFile($qcfile, $name);
    }

    if ((-e $errfile) && (! -z $errfile)) {
        my @n       = split "/", $errfile;
        my $l       = scalar(@n);
        my $name    = "$n[$l-3]-$n[$l-2]-error.txt";

        attachFile($errfile, $name);
    }
}


########################################
#
#  All done.  DO NOT add another boundary - this results in a bogus empty attachment.
#
#print "--Bri_Says_This_Is_The_Boundary\n";








sub attachFile ($$) {
    my $fileName   = shift @_;  #  Name of file on disk
    my $attachName = shift @_;  #  Name of file in email

    return  if (! -e $fileName);
    return  if (-z   $fileName);

    print "--Bri_Says_This_Is_The_Boundary\n";
    print "Content-Type: application/octet-stream; name=\"$attachName\"\n";
    print "Content-Transfer-Encoding: base64\n";
    print "Content-Disposition: attachment; filename=\"$attachName\"\n";
    print "\n";

    my @message;
    my $message;

    #  Winblows wants \r\n, and since that's the primary consumer of these emails, we put it in.

    open(F, "< $fileName");
    while (<F>) {
        chomp;
        push @message, "$_\r\n";
    }
    close(F);

    $message = join "", @message;

    local($/) = undef;
    print encode_base64($message);
}
