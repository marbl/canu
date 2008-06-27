#!/usr/bin/perl -w
use strict;
use AMOS::AmosLib;

my $USAGE = "stripCGBReads.pl list.txt < cgb.orig > cgb\n";


my $list = shift @ARGV or die $USAGE;
open LIST, "< $list" or die "Can't open $list ($!)\n";

my %skipids;
while (<LIST>)
{
  chomp;
  $skipids{$_} = -1;

  print STDERR "Marked $_ to strip\n";
}


while (my $rec = getRecord(\*STDIN))
{
  my ($id, $fields, $recs) = parseRecord($rec);

  if ($id eq "IUM")
  {
    my $unitignum = $fields->{acc};
    foreach my $frg (@$recs)
    {
      my ($fid, $ffields, $rrecs) = parseRecord($frg);

      if ($fid eq "IMP")
      {
        my $mid = $ffields->{mid};
        my $con = $ffields->{con};

        if (exists $skipids{$mid})
        {
          my $name = $ffields->{src} || "UNDEF";
          chomp $name;

          print STDERR "Stripping read $mid ($name) occurring in unitig $unitignum contained by $con\n";
          $skipids{$mid} = $con;

          $fields->{nfr}--;

          #$rec = sprintf("{IUM\nacc:%d\nsrc:\n%s.\ncov:%s\nsta:%s\nabp:%d\nbbp:%d\nlen:%d\ncns:\n.\nqlt:\n.\nfor:%d\nnfr:%d\n",
          #               $fields->{acc}, $fields->{src}, $fields->{cov}, $fields->{sta}, $fields->{abp}, $fields->{bbp}, $fields->{len},$fields->{for},$fields->{nfr});

          $rec = sprintf("{IUM\nacc:%d\ncov:%s\nsta:%s\nabp:%d\nbbp:%d\nlen:%d\ncns:\n.\nqlt:\n.\nfor:%d\nnfr:%d\n",
                         $fields->{acc}, $fields->{cov}, $fields->{sta}, $fields->{abp}, $fields->{bbp}, $fields->{len},$fields->{for},$fields->{nfr});


          foreach my $f (@$recs)
          {
            my ($fid, $ff, $rrecs) = parseRecord($f);
            next if $ff->{mid} eq $mid;

            if ($ff->{con} eq $mid)
            {
              print STDERR "Read $ff->{mid} was contained by $mid, repointing\n";
              $ff->{con} = $skipids{$mid};
            }

            #$rec .= sprintf("{IMP\ntyp:%s\nmid:%s\ncon:%s\nsrc:\n%s.\npos:%s\ndln:%s\ndel:\n}\n",
            #                $ffields->{typ}, $ff->{mid}, $ff->{con}, $ff->{src}, $ff->{pos}, $ff->{dln}, $ff->{del});

            $rec .= sprintf("{IMP\ntyp:%s\nmid:%s\ncon:%s\npos:%s\ndln:%s\ndel:\n}\n",
                            $ffields->{typ}, $ff->{mid}, $ff->{con}, $ff->{pos}, $ff->{dln}, $ff->{del});

          }

          $rec .= "}\n";
        }
      }
    }
  }

  print $rec;
}
