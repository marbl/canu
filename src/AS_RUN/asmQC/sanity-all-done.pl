#!/usr/bin/perl

use MIME::QuotedPrint;
use MIME::Base64;
use HTML::Entities;
use strict;

if (scalar(@ARGV) != 5) {
    die "wrong args.\n";
}

my $bindir    = shift @ARGV;
my $wrkdir    = shift @ARGV;
my $thisdate  = shift @ARGV;
my $lastdate  = shift @ARGV;
my $addresses = shift @ARGV;
my $timenow   = localtime();


print "Subject: CAtest $thisdate\n";
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
print "<FONT SIZE=+2>" . encode_entities("Results for $thisdate (Finished at $timenow).") . "</FONT>\n";
print "</P>\n";
print "\n";

print "<P>\n";

open(R, "cat $wrkdir/$thisdate/*/sanity-result.out |");
while (<R>) {
    chomp;
    s!SUCCESS!<FONT COLOR=\"0x00ff00\">SUCCESS</FONT>!;
    s!FAILURE!<FONT COLOR=\"0xff0000\">FAILURE</FONT>!;
    print "$_<BR>\n";
}
close(R);

print "</P>\n";

print "<HR>\n";
print "<P>\n";
print "<FONT SIZE=+2>Changes for kmer.</FONT>\n";
print "</P>\n";
open(F, "< $wrkdir/$thisdate/wgs/kmer.updates") or die;
while (<F>) {
    chomp;
    print encode_entities($_) . "<BR>\n";
}
close(F);

print "<HR>\n";
print "<P>\n";
print "<FONT SIZE=+2>Changes for wgs-assembler.</FONT>\n";
print "</P>\n";
open(F, "< $wrkdir/$thisdate/wgs/src.updates") or die;
while (<F>) {
    chomp;
    print encode_entities($_) . "<BR>\n";
}
close(F);

print "<HR>\n";
print "<P>\n";
print "<FONT SIZE=+2>Build Results for kmer.</FONT>\n";
print "</P>\n";
open(F, "< $wrkdir/$thisdate/wgs/kmer/make.err") or die;
while (<F>) {
    chomp;
    print encode_entities($_) . "<BR>\n";
}
close(F);

print "<HR>\n";
print "<P>\n";
print "<FONT SIZE=+2>Build Results for wgs-assembler.</FONT>\n";
print "</P>\n";
open(F, "< $wrkdir/$thisdate/wgs/src/make.err") or die;
while (<F>) {
    chomp;
    print encode_entities($_) . "<BR>\n";
}
close(F);

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
#  Dump the QC attachments
#
open(R, "ls $wrkdir/$thisdate/*/sanity-error.out |");
while (<R>) {
    chomp;

    my $errfile = $_;
    my $qcfile  = $_;

    $qcfile =~ s/sanity-error.out/sanity-qc.out/;
    
    #  Winblows wants \r\n, and since that's the primary consumer of these emails, we
    #  put it in.

    if (! -z $qcfile) {
        my @n       = split "/", $qcfile;
        my $l       = scalar(@n);
        my $name    = "$n[$l-3]-$n[$l-2]-qc.txt";

        print "--Bri_Says_This_Is_The_Boundary\n";
        print "Content-Type: application/octet-stream; name=\"$name\"\n";
        print "Content-Transfer-Encoding: base64\n";
        print "Content-Disposition: attachment; filename=\"$name\"\n";
        print "\n";

        my @message;
        my $message;

        open(F, "< $qcfile") or die;
        while (<F>) {
            chomp;
            push @message, "$_\r\n";
        }
        close(F);

        $message = join "", @message;

        local($/) = undef;
        print encode_base64($message);
    }

    if (! -z $errfile) {
        my @n       = split "/", $errfile;
        my $l       = scalar(@n);
        my $name    = "$n[$l-3]-$n[$l-2]-error.txt";

        print "--Bri_Says_This_Is_The_Boundary\n";
        print "Content-Type: application/octet-stream; name=\"$name\"\n";
        print "Content-Transfer-Encoding: base64\n";
        print "Content-Disposition: attachment; filename=\"$name\"\n";
        print "\n";

        my @message;
        my $message;

        open(F, "< $errfile") or die;
        while (<F>) {
            chomp;
            push @message, "$_\r\n";
        }
        close(F);

        $message = join "", @message;

        local($/) = undef;
        print encode_base64($message);
    }

}


########################################
#
#  All done.  DO NOT add another boundary - this results in a bogus empty attachment.
#
#print "--Bri_Says_This_Is_The_Boundary\n";
