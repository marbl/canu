#!/usr/local/bin/perl
use IO::File;
use File::Basename;
use TIGR::Foundation;
use TIGR::AsmLib;
use strict;

my $MY_VERSION = " Version 1.1 (Build " . (qw/$Revision: 1.2 $/ )[1] . ")";
my @DEPENDS =
(
  "TIGR::Foundation"
);



my $MY_HELPTEXT = qq~
    ca2scaff -i file.asm -o prefix

    Create gaps, oo, sum, and scaff files of the Celera Assembler scaffold
    within the asm file.

    .scaff format:

      >scaffid numcontigs scaffbases scaffspan
      contigid orientation contiglen gapsize gapstdeviation
      ...

    .gaps format:

      >scaffid
      contigid contigid gapsize gapstdeviation
      ...

    .oo format:

      >scaffid
      contigid orientation
      ...

    .sum format

      scaffid numcontigs scaffbases scaffspan
      ...

    Options
    -------
    -s Only print the scaff file

    -i <asmfile> Specfiy input asm file
    -o <prefix>  Specify prefix to output files
~;

my $base = new TIGR::Foundation;

if (! defined $base){
    print STDERR "Nasty error, hide!\n";
    exit(1);
}


$base->setHelpInfo($MY_HELPTEXT);
$base->addDependInfo(@DEPENDS);
$base->setVersionInfo($MY_VERSION);

my $infile;
my $outfile;
my $scaffonly = 0;

my $err = $base->TIGR_GetOptions("i=s" => \$infile,
                                 "o=s" => \$outfile,
                                 "s"   => \$scaffonly);

$base->bail("Command line parsing failed.  See -h option")
  if ($err == 0);

if (! defined $infile)
{
  if ($#ARGV < 0)
  {
    $base->bail("Must specify an input file name.  See -h option");
  } else {
    $infile = $ARGV[0];
  }
}

$base->bail("Must specify an output prefix.  See -h option")
  if (! defined $outfile);

my $oofile    = "$outfile.oo";
my $sumfile   = "$outfile.sum";
my $gapfile   = "$outfile.gaps";
my $scafffile = "$outfile.scaff";

open(SCAFF, ">$scafffile") || $base->bail("Cannot open $scafffile: $!\n");

if (!$scaffonly)
{
  open(OO,    ">$oofile")    || $base->bail("Cannot open $oofile: $!\n");
  open(SUM,   ">$sumfile")   || $base->bail("Cannot open $sumfile: $!\n");
  open(GAP,   ">$gapfile")   || $base->bail("Cannot open $gapfile: $!\n");
}

open(IN, $infile) || $base->bail("Cannot open $infile: $!");


my %ctglen;
my %ctglinks;
my $record;

while ($record = getCARecord(\*IN))
{
  my ($type, $fields, $recs) = parseCARecord($record);

  my $id = getCAId($$fields{acc});

  if ($type eq "CCO")
  {
    my $contiglen = $$fields{"len"};
    while ($$fields{"cns"} =~ /-/g){ $contiglen--; }

    $ctglen{$id} = $contiglen;
  }
  elsif ($type eq "SCF")
  {
    my $scflen = 0;
    my $scfspan = 0;
    my $nctg = $$fields{"noc"} + 1;

    ## Do a quick pass to compute scaffold len,span
    for (my $i = 0; $i <= $#$recs; $i++)
    {
      my ($sid, $sfs, $srecs) = parseCARecord($$recs[$i]);

      if ($sid eq "CTP")
      {
        my $c1 = $$sfs{"ct1"};
        my $c2 = $$sfs{"ct2"};

        if ($i == 0)
        {
          $scflen  += $ctglen{$c1};
          $scfspan += $ctglen{$c1};
        }

        if ($nctg != 1)
        {
          my $scaffgap = $$sfs{mea};

          $scflen  += $ctglen{$c2};
          $scfspan += $ctglen{$c2} + $scaffgap;
        }
      }
    }

    $scfspan = int($scfspan);

    ## Now print the scaffold and contigs
    print SCAFF ">$id $nctg $scflen $scfspan\n";

    if (!$scaffonly)
    {
      print OO    ">$id\n";
      print GAP   ">$id\n";
    }

    for (my $i = 0; $i <= $#$recs; $i++)
    {
      my ($sid, $sfs, $srecs) = parseCARecord($$recs[$i]);

      if ($sid eq "CTP")
      {
        my $c1  = $$sfs{"ct1"};
        my $c2  = $$sfs{"ct2"};
        my $ori = $$sfs{"ori"};

        my ($o1, $o2);
        if ($ori eq "N" || $ori eq "I"){ $o1 = "BE"; }
        else                           { $o1 = "EB"; }

        if ($ori eq "N" || $ori eq "O"){ $o2 = "BE"; }
        else                           { $o2 = "EB"; }

        if ($i == 0)
        {
          if (!$scaffonly)
          {
            print OO "$c1 $o1\n";
          }
        }

        if ($nctg != 1)
        {
          if (!$scaffonly)
          {
            print OO "$c2 $o2\n";
          }

          my $scaffgap = $$sfs{mea};
          my $scaffstd = $$sfs{std};

          print SCAFF "$c1 $o1 $ctglen{$c1} $scaffgap $scaffstd\n";

          if ($i == $#$recs)
          {
            print SCAFF "$c2 $o2 $ctglen{$c2} 0 0.0\n";
          }

          if (!$scaffonly)
          {
            print GAP "$c1 $c2 $scaffgap $scaffstd\n";
          }
        }
        else
        {
          print SCAFF "$c1 $o1 $ctglen{$c1} 0 0.0\n";
        }
      } # if CTP
    } # for each rec

    if (!$scaffonly)
    {
      print SUM "$id $nctg $scflen ", int($scfspan), "\n";
    }
  } # if $type = SCF
}

close(SCAFF);
close(IN);

if (!$scaffonly)
{
  close(SUM);
  close(GAP);
  close(OO);
}

exit(0);
