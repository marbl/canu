#!/usr/bin/perl

#   Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
#   all rights reserved. Authored by: Sergey Koren
#   

#   This Software was prepared for the Department of Homeland Security
#   (DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
#   part of contract HSHQDC-07-C-00020 to manage and operate the National
#   Biodefense Analysis and Countermeasures Center (NBACC), a Federally
#   Funded Research and Development Center.
#   
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#   
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   
#   * Neither the name of the Battelle National Biodefense Institute nor
#     the names of its contributors may be used to endorse or promote
#     products derived from this software without specific prior written
#     permission.
#   
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
###########################################################################
#
#  Read in fragments from fastq-format sequence and quality files,
#  correct to pacbio fragments.
#

use strict;

use Config;  #  for @signame
use FindBin;
use Cwd;
use Carp;
use FileHandle;

use POSIX qw(ceil floor sys_wait_h);

my %global;

sub getGlobal ($) {
    my $var = shift @_;
    if (!exists($global{$var})) {
        return undef;
    }
    return($global{$var});
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;

    #  If no value, set the field to undefined, the default for many of the options.

    $val = undef  if ($val eq "");

    $global{$var} = $val;
    $val =~ s/\\\"/\"/g;
    $val =~ s/\"/\\\"/g;
    $val =~ s/\\\$/\$/g;
    $val =~ s/\$/\\\$/g;
}

sub readFasta($$$) {
   my $file = shift @_;
   my $prefix = shift @_;
   my $hashref = shift @_;

   my $fastaSeq = "";
   my $id = "";
   open(CNS, "< $file") or die("Couldn't open '$file'", undef);
   while (<CNS>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);

      if (m/\>/) {
         if ($fastaSeq ne "") {
            $fastaSeq =~ s/^\s+|\s+$//g;
            $hashref->{$id} = $fastaSeq;
         }
         $fastaSeq = "";
         $id = (split '\s+')[0];
         $id =~ s/>//g;
         if ($id =~ m/utg_/) {
            $id = (split '/', $id)[0];
            $id =~ s/utg_//g;
         } else {
            $id = (split ',', $id)[1];
         }
         $id = "$prefix.$id";
      } else {
         $fastaSeq = $fastaSeq . $_;
      }
   }
   if ($fastaSeq ne "") {
      $fastaSeq =~ s/^\s+|\s+$//g;
      $hashref->{$id} = $fastaSeq;
   }
   close(CNS);
}

sub process($$$$) {
   my $cns = shift @_;
   my $raw = shift @_;
   my $lay = shift @_;
   my $output = shift @_;
   my %seqs = {};

   readFasta($cns, "utg", \%seqs);
   readFasta($raw, "seq", \%seqs);
   my $count = (scalar keys %seqs);
   print "Loaded $count seqs\n";

   my $qltChar = chr(getGlobal("quality") + ord('0'));
   # we scan file twice, could probably be done once 
   # first scan figure out which sequences have only one read and get their cns from the read or full cns
   my %cnsHash = {};
   my $id = "";
   open(LAY, "<$lay") or die ("Couldn't open '$lay'", undef);
   while (<LAY>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);
      my $type = (split '\s+')[0];
      my $val = (split '\s+')[1];
      if ($type eq "unitig") {
         $id = $val;
      } elsif ($type eq "contig") {
         $id = $val;
         $cnsHash{$id} = "utg.$id";
      } elsif ($type eq "data.num_frags" && $val > 1) {
         $cnsHash{$id} = "utg.$id";
      } elsif ($type eq "FRG" && !defined($cnsHash{$id})) { # this is a single frg utg
         $cnsHash{$id} = "seq." . (split '\s+')[4];
      }
   }
   close(LAY);

   my $fastaSeq = "";
   my $utgID = "";
   my $isUTG = 1;
   open(OUT, ">$output") or die ("Couldn't open '$output'", undef);
   open(LAY, "<$lay") or die ("Couldn't open '$lay'", undef);
   while (<LAY>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);
      my $type = (split '\s+')[0];
      my $val = (split '\s+')[1];
      if ($type eq "unitig") {
         # close previous unitig
         if ($fastaSeq ne "") {
            print OUT "UTG type X ident\t$utgID position\t0\t" . length($fastaSeq) . " num_instances 0\n";
         }
         my $cnsId = $cnsHash{$val};
         # if unknown unitig encountered due to it being too short delete it from the tigStore
         if (!defined($seqs{$cnsId})) {
            print STDERR "Unknown consensus for unitig $val\n";
            $fastaSeq="";
            print OUT "$_\n";
	    print OUT "len 0\n";
	    print OUT "cns\n";
	    print OUT "qlt\n";
	    print OUT "data.unitig_coverage_stat\t1.000000\n";
	    print OUT "data.unitig_microhet_prob\t1.000000\n";
	    print OUT "data.unitig_status\tX\n";
	    print OUT "data.unitig_suggest_repeat\tF\n";
	    print OUT "data.unitig_suggest_unique\tF\n";
	    print OUT "data.unitig_force_repeat\tF\n";
	    print OUT "data.unitig_force_unique\tF\n";
            print OUT "data.contig_status\tU\n";
            print OUT "data.num_frags\t0\n";
            print OUT "data.num_unitigs\t0\n";
         } else {
            $fastaSeq = uc $seqs{$cnsId};
            $utgID = $val;
            print OUT "$_\n";
         }
      } elsif ($type eq "contig") {
         # close previous contig
         my $cnsId = $cnsHash{$val};
         die "Unknown consensus for contig $val\n" if (!defined($seqs{$cnsId}));
         $fastaSeq = uc $seqs{$cnsId};
         $utgID = $val;
         $isUTG = 0;
         setGlobal("type", "-cp");
         print OUT "$_\n";
      } elsif ($fastaSeq ne "") {
         if ($type eq "len") {
            print OUT "len " . length($fastaSeq) . "\n";
         } elsif ($type eq "cns") {
            print OUT "cns $fastaSeq\n";
         } elsif ($type eq "qlt") {
            my $qlt = $fastaSeq;
            $qlt =~ s/./$qltChar/g;
            print OUT "qlt $qlt\n";
         } elsif ($type eq "data.num_unitigs" && $val < 1) {
            print OUT "data.num_unitigs\t1\n";
         } else {
            if ($type eq "FRG") {
               my @tokenized =  split '\s+';
               my $min = ($tokenized[13] < $tokenized[14] ? $tokenized[13] : $tokenized[14]);
               my $max = ($tokenized[13] > $tokenized[14] ? $tokenized[13] : $tokenized[14]);
               my $fwd = ($tokenized[13] < $tokenized[14] ? "+" : "-");

               if ($fastaSeq ne "" && $max > length($fastaSeq)) {
                  $max = length($fastaSeq);
                  if ($fwd eq "+") { 
                     $tokenized[14] = $max;
                  } else {
                     $tokenized[13] = $max;
                  }
               }
               if ($fastaSeq ne "" && $min > length($fastaSeq)) {
                  $min = length($fastaSeq);
                  if ($fwd eq "+") {
                     $tokenized[13] = $min;
                  } else {
                     $tokenized[14] = $min;
                  }
               }
               $_ = "";
               for (my $i = 0; $i < scalar @tokenized; $i++) {
                  $_ = $_ . "\t" . $tokenized[$i];
               }
               
            }
            print OUT "$_\n";
         }
      }
   }
   if ($fastaSeq ne "" && $isUTG) {
      print OUT "UTG type X ident\t$utgID position\t0\t" . length($fastaSeq) . " num_instances 0\n";
   }

   close(LAY);
   close(OUT);
}

my $err = 0;
setGlobal("prefix", "tmp");
setGlobal("quality", 60);
setGlobal("type", "-up");
my $input = undef;
my $seq = undef;
my $output = undef;
my $lay = undef;
my $partition = undef;
my $version = 1;

# finally, command-line parameters take precedence
while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if ($arg eq "-prefix") {
        setGlobal("prefix", shift @ARGV);

    } elsif ($arg eq "-partition") {
        $partition = shift @ARGV;

    } elsif ($arg eq "-path") {
        setGlobal("path", shift @ARGV);

    } elsif ($arg eq "-input") {
       $input = shift @ARGV;

    } elsif ($arg eq "-lay") {
       $lay = shift @ARGV;

    } elsif ($arg eq "-output") {
        $output = shift @ARGV;

    } elsif ($arg eq "-sequence") {
        $seq = shift @ARGV;

    } elsif ($arg eq "-quality") {
        setGlobal("quality", shift @ARGV);

    } elsif ($arg eq "-version") {
        $version = shift @ARGV;

    } else {
       print STDERR "Unknown parameter " + $err + "\n";
       $err++;
    }
 }

if ($err > 0 || !defined($input) || !defined($output) || !defined($lay)) {
   print STDERR "Error: no valid input cns specified\n" if (!defined($input));
   print STDERR "Error: no valid output fasta specified\n" if (!defined($output));
   print STDERR "Error: no valid input lay specified\n", if (!defined($lay));
   print STDERR "Invalid arguments specified\n";

   print STDERR "usage: $0 [options] -input <input cns> -output <output layout> -lay <input layout>\n";
   print STDERR "  -path                 Path to smrtportal installation. Required if smrtportal is not in path.\n";
   print STDERR "  -prefix               Prefix for stores to write to\n";
   exit(1);
}
# find tigStore
my $TIG = getGlobal("path");
if (!defined($TIG) || $TIG eq "") {
   print STDERR "Error, could not find tigStore in path. Please update your path or specify a location using -path\n";
   exit(1);
}

$output="$output.$partition.cnslay";
# clear output file
open(O, ">", $output);
close(O);

process($input, $seq, $lay, $output);
# update store
my $prefix = getGlobal("prefix");
my $type = getGlobal("type");
system("$TIG/tigStore -g $prefix.gkpStore -t $prefix.tigStore $version $type $partition -N -A -R $output");
if ($? != 0) {
   die "Error: tigStore could not run successfully"
}
# cleanup
system("rm -f $output");
