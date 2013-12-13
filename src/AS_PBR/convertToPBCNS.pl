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


sub processLayout($$$$) {
   my $output = shift @_;
   my $id = shift @_;
   my $cns = shift @_;
   my $hashref = shift @_;
   my %hash = %$hashref;
   my $prefix = getGlobal("prefix");
   my $path = getGlobal("smrtpath");
   my $threads = getGlobal("threads");
   my $cov = getGlobal("coverage");
   my $length = getGlobal("length");

   open(CNS, "> $prefix.cns.fasta") or die("Couldn't open '$prefix.cns.fasta'", undef);
   $cns =~ s/^\s+|\s+$//g;
   print CNS ">$id\n";
   print CNS "$cns\n";
   close(CNS);

   open(ALN, "> $prefix.aln.fasta") or die ("Couldn't open '$prefix.aln.fasta'", undef);
   for (keys %hash) {
      print ALN ">$_\n";
      my $aln = $hash{$_};
      $aln =~ s/^\s+|\s+$//g;
      print ALN "$aln\n";
   }
   close(ALN);

   system("$path/blasr $prefix.aln.fasta $prefix.cns.fasta -m 5 -nproc $threads > $prefix.m5 2>/dev/null");
   if ($? != 0) { 
      die "Error: blasr could not run successfully"
   }
   system("$path/pbdagcon -v -t 0 -m $length -j $threads -c $cov $prefix.m5 >> $output 2>/dev/null");
   if ($? != 0) {
      die "Error: pbdagon could not run successfully"
   }
}
 
sub processLayouts($$) {
   my $input = shift @_;
   my $output = shift @_;

   my $lastEnd = 0;
   my $cns = "";
   my $lastID = "";
   my %toAlign = {};
   open(F, "< $input") or die("Couldn't open '$input'", undef);
   while (<F>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);
      my @tokenized =  split '\s+';
      my $hashLen = scalar keys %toAlign;

      if ($lastID eq "" || $lastID ne $tokenized[1]) {
         if ($lastID ne "") { processLayout($output, $lastID, $cns, \%toAlign); } 
         for (keys %toAlign) { delete $toAlign{$_}; }
         $cns = ' ' x ($tokenized[3] + 1);
         $lastEnd = 0;
      }
      my $end = $tokenized[5];
      my $start = $tokenized[4];
      if ($end >= $lastEnd) {
         if ($end - $start > length($tokenized[7])) {
            $end = $start + length($tokenized[7]);
         }
         $lastEnd = $end;
         substr($cns, $start, $end-$start) = $tokenized[7];
      }
      $toAlign{$tokenized[0]} = $tokenized[6];
      $lastID = $tokenized[1];
      
   }
   if ($lastID ne "") {
      processLayout($output, $lastID, $cns, \%toAlign);
   }
   close(F);
}

my $err = 0;
setGlobal("prefix", "tmp");
setGlobal("coverage", 0);
setGlobal("threads", 1);
my $input = undef;
my $output = undef;

# finally, command-line parameters take precedence
while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;

    if ($arg eq "-prefix") {
        setGlobal("prefix", shift @ARGV);

    } elsif ($arg eq "-coverage") {
        my $cov = shift @ARGV;
        if ($cov < 0) { $cov = 0; }
        setGlobal("coverage", $cov);

    } elsif ($arg eq "-threads") {
        my $threads = shift @ARGV;
        if ($threads < 1) { $threads = 1; }
        setGlobal("threads", $threads);
 
    } elsif ($arg eq "-length") {
        my $length = shift @ARGV;
        if ($length < 1) { $length = 1; }
        setGlobal("length", $length);

    } elsif ($arg eq "-path") {
        setGlobal("smrtpath", shift @ARGV);

    } elsif ($arg eq "-input") {
       $input = shift @ARGV;

    } elsif ($arg eq "-output") {
        $output = shift @ARGV;

    } else {
       print STDERR "Unknown parameter " + $err + "\n";
       $err++;
    }
 }

if ($err > 0 || !defined($input) || !defined($output)) {
   print STDERR "Error: no valid input lay specified\n" if (!defined($input));
   print STDERR "Error: no valid output fasta specified\n" if (!defined($output));
   print STDERR "Invalid arguments specified\n";

   print STDERR "usage: $0 [options] -input <input layout> -output <output fasta>\n";
   print STDERR "  -coverage              Minimum coverage in a consensus region to keep, default 0.\n";
   print STDERR "  -threads              Number of threads to use for generating consensus, default 1.\n";
   print STDERR "  -path                 Path to smrtportal installation. Required if smrtportal is not in path.\n";
   print STDERR "  -prefix               Prefix for temporary files created by this job\n";
   print STDERR "  -length               Minimum sequence length to keep\n";
   exit(1);
}

# find blasr
my $BLASR = getGlobal("smrtpath");
if (!defined($BLASR) || $BLASR eq "") {
   # try to use path
   my $blasrPath = `which blasr`;
   chomp $blasrPath;
   my @t = split '/', "$blasrPath";
   pop @t;                      #  blasr 
   $BLASR = join '/', @t;  #  path to the assembler
   setGlobal("smrtpath", $BLASR);
}
# if we really can't find it just give up
if (! -e "$BLASR/blasr" || ! -e "$BLASR/pbdagcon") {
   print STDERR "Error, could not find blasr/smrtportal in path. Please update your path or specify a location using -path\n";
   exit(1);
}

# clear output file
open(O, ">", $output);
close(O);
processLayouts($input, $output);
# cleanup
system("rm -f " . getGlobal("prefix") . "*");
