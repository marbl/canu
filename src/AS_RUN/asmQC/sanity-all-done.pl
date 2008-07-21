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
my $lastdate  = shift @ARGV;
my $thisdate  = shift @ARGV;
my $addresses = shift @ARGV;
my $timenow   = localtime();


open(F, "ls $wrkdir/$thisdate/*/assembly-done.out |");
my @asmResults = <F>;
chomp @asmResults;
close(F);


print "Subject: CAtest $thisdate\n";
print "To:      $addresses\n";
print "Content-Type: multipart/mixed; boundary=Bri_Says_This_Is_The_Boundary\n";  # alternative if no atachments
print "\n";


########################################
#
#  Dump the ASCII part
#

if (0) {
print "--Bri_Says_This_Is_The_Boundary\n";
print "Content-Type: text/plain; charset=\"iso-8859-1\"\n";
print "Content-Transfer-Encoding: quoted-printable\n";
print "\n";

print "Results for $thisdate (Finished at $timenow).\n";
print "\n";

foreach my $asmres (@asmResults) {
    open(F, "< $asmres") or die;
    while (<F>) {
        chomp;
        if (m/Assembly\sresult/) {
            print "$_\n";
        }
    }
    close(F);
}
}


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
foreach my $asmres (@asmResults) {
    open(F, "< $asmres") or die;
    while (<F>) {
        chomp;
        if (m/Assembly\sresult/) {
            s!SUCCESS!<FONT COLOR=\"0x00ff00\">SUCCESS</FONT>!;
            s!FAILURE!<FONT COLOR=\"0xff0000\">FAILURE</FONT>!;
            print "$_<BR>\n";
        }
    }
    close(F);
}
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
#  Dump the QC attachments
#

foreach my $asmres (@asmResults) {
    my @n    = split "/", $asmres;
    my $l    = scalar(@n);
    my $name = "$n[$l-3]-$n[$l-2]-qc.txt";

    print "--Bri_Says_This_Is_The_Boundary\n";
    print "Content-Type: application/octet-stream; name=\"$name\"\n";
    print "Content-Transfer-Encoding: base64\n";
    print "Content-Disposition: attachment; filename=\"$name\"\n";
    print "\n";

    #  Winblows wants \r\n, and since that's the primary consumer of these emails, we
    #  put it in.

    my @message;
    my $message;

    open(F, "< $asmres") or die;
    while (<F>) {
        chomp;
        push @message, "$_\r\n";
    }
    close(F);

    $message = join "", @message;

    local($/) = undef;
    print encode_base64($message);
}


########################################
#
#  All done.
#

#print "--Bri_Says_This_Is_The_Boundary\n";


__END__




foreach my $asmres (@asmResults) {
    open(F, "< $asmres") or die;

    $_ = <F>;
    chomp;

    print "<HR>\n";
    print "<HR>\n";
    print "<P>\n";
    print "<FONT SIZE=+2>" . encode_entities($_) . "</FONT>\n";
    print "</P>\n";

    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;

    print "<TABLE BORDER=1>\n";
    print "<TR><TH></TH><TH>ref</TH><TH>$lastdate</TH><TH>$thisdate</TH>\n";
    while (<F>) {
        my @v = split '\s+', $_;
        print "<TR><TD>" . encode_entities($v[0]) . "</TD><TD>" . encode_entities($v[1]) . "</TD><TD>" . encode_entities($v[2]) . "</TD><TD>" . encode_entities($v[2]) . "</TD></TR>\n";
    }
    print "</TABLE>\n";
}
