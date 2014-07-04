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
my $layCounter = 0;

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

    if ($var eq "consensus") {
       if ($val eq "pbdagcon") {
          setGlobal("batch", 1000);
          setGlobal("runBlasr", 1);
       } elsif ($val eq "falcon" ){
          setGlobal("batch", 1000);
          setGlobal("runBlasr", 0);
       } elsif ($val eq "fast") {
          setGlobal("batch", 1);
          setGlobal("runBlasr", 0);
       } elsif ($val eq "pbutgcns") {
          setGlobal("batch", 1);
          setGlobal("runBlasr", 0);
       } else {
          die "Unknown consensus type $val specified.\n";
       }
    }
}

sub flushFalcon($) {
   my $output = shift @_;
   my $consensus = getGlobal("consensus");

   if ($consensus ne "falcon") { return; }

   my $threads = getGlobal("threads");
   my $cov = getGlobal("coverage");
   my $path = getGlobal("smrtpath");
   my $python = getGlobal("falconpath");
   my $prefix = getGlobal("prefix");
   my $fileName = "$prefix.cns.in";

   open(CNS, ">> $fileName") or die("Couldn't open '$fileName'", undef);
   print CNS "- -\n";
   close(CNS);

   my $libPath = "$python";
   $libPath =~ s/bin/lib/g;
   my $smrtPortalPath = "$path";
   $smrtPortalPath =~ s/bin/lib/g;
   my $pythonPath = "$path:$libPath/python2.7/site-packages:$smrtPortalPath/python2.7:\$PYTHONPATH";
   system("export PYTHONPATH=$pythonPath && cat $fileName | $python/falcon_sense.py --trim 1 --local_match_count_threshold 3 --min_cov $cov --n_core $threads >> $output 2> $prefix.err");

   if ($? != 0) {
      die "Error: falcon_sense could not run successfully"
   }
}

sub flushPBDAGCON($) {
   my $output = shift @_;
   my $consensus = getGlobal("consensus");

   if ($consensus ne "pbdagcon") { return; }

   my $threads = getGlobal("threads");
   my $cov = getGlobal("coverage");
   my $path = getGlobal("smrtpath");
   my $length = getGlobal("length");
   my $python = getGlobal("falconpath");
   my $prefix = getGlobal("prefix");

   if (! -e "$prefix.m5" || ! -s "$prefix.m5") { return; }
   system("$path/pbdagcon -t 0 -m $length -j $threads -c $cov $prefix.m5 >> $output 2>/dev/null");
   if ($? != 0) {
      die "Error: pbdagon could not run successfully"
   }
}

sub flushPBUTGCNS($) {
   my $output = shift @_;
   my $consensus = getGlobal("consensus");

   if ($consensus ne "pbutgcns") { return; }

   my $threads = getGlobal("threads");
   my $cov = getGlobal("coverage");
   my $path = getGlobal("smrtpath");
   my $length = getGlobal("length");
   my $python = getGlobal("falconpath");
   my $prefix = getGlobal("prefix");

   if (! -e "$prefix.cns.in" || ! -s "$prefix.cns.in") { return; }
   system("$path/pbutgcns -j $threads $prefix.cns.in >> $output 2>/dev/null");
   if ($? != 0) {
      die "Error: pbutgcns could not run successfully"
   }
}

sub processLayout($$$$$) {
   my $output = shift @_;
   my $id = shift @_;
   my $cns = shift @_;
   my $hashref = shift @_;
   my $sensitive = shift @_;
   my %hash = %$hashref;
   my $prefix = getGlobal("prefix");
   my $path = getGlobal("smrtpath");
   my $threads = getGlobal("threads");
   my $cov = getGlobal("coverage");
   my $length = getGlobal("length");
   my $consensus = getGlobal("consensus");

   my $alnCount = scalar keys %hash;
   $cns =~ s/^\s+|\s+$//g;

   # just output our frakenconsensus
   if ($consensus eq "fast") {
      open(OUT, ">> $output") or die("Couldn't open '$output'", undef);
      print OUT ">$id\n";
      print OUT "$cns\n";
      close(OUT);
      return;
   }

   #  now we format our data for whatever consensus was selected
   if ($alnCount == 0) {
      return;
   } elsif ($alnCount == 1) {
      if ($cov <= 1) {
         open(OUT, ">> $output") or die("Couldn't open '$output'", undef);
         print OUT ">$id\n";
         print OUT "$cns\n";
         close(OUT);
      } else {
         # skip low coverage stuff
         return;
       }
   } elsif (length($cns) < $length) {
      # skip short stuff
      return;
   } elsif ($consensus eq "pbdagcon") {
         open(CNS, "> $prefix.cns.fasta") or die("Couldn't open '$prefix.cns.fasta'", undef);
         print CNS ">$id\n";
         print CNS "$cns\n";
         close(CNS);

         open(ALN, "> $prefix.aln.fasta") or die ("Couldn't open '$prefix.aln.fasta'", undef);
         for (keys %hash) {
            print ALN ">$_\n";
            my $str = $hash{$_};
            my @tokenized = split('\s+', $str);
            my $aln =  $tokenized[2];
            $aln =~ s/^\s+|\s+$//g;
            print ALN "$aln\n";
         }
         close(ALN);

         if ($sensitive) {
            system("$path/blasr $prefix.aln.fasta $prefix.cns.fasta -sdpDel 5 -minMatch 10 -maxMatch 16 -bestn 1 -m 5 -nproc $threads >> $prefix.m5 2>/dev/null");
         } else {
            system("$path/blasr $prefix.aln.fasta $prefix.cns.fasta -maxMatch 16 -bestn 1 -m 5 -nproc $threads >> $prefix.m5 2>/dev/null");
         }
         if ($? != 0) { 
            die "Error: blasr could not run successfully"
         }
         $layCounter++;
         if ($layCounter > getGlobal("batch")) {
            flushPBDAGCON($output);
            open(CNS, "> $prefix.m5") or die("Couldn't open '$prefix.m5'", undef);
            close(CNS);
            $layCounter = 0;
         }
   } elsif ($consensus eq "pbutgcns") {
      if ($layCounter == 0) {
         open(CNS, "> $prefix.cns.in") or die ("Couldn't open '$prefix.cns.in'", undef);
      } else {
         open(CNS, ">> $prefix.cns.in") or die ("Couldn't open '$prefix.cns.in'", undef);
      }
      print CNS "$id $cns\n";
      for (keys %hash) {
         my $str = $hash{$_};
         my @tokenized = split('\s+', $str);
         my $aln = $tokenized[2];
         $aln =~ s/^\s+|\s+$//g;
         print CNS "$_ $tokenized[0] $tokenized[1] $aln\n";
      }
      close(CNS);
      $layCounter++;
      if ($layCounter >= getGlobal("batch")) {
         flushPBUTGCNS($output);
         open(CNS, "> $prefix.cns.in") or die("Couldn't open '$prefix.cns.in'", undef);
         close(CNS);
         $layCounter = 0;
      }
   } elsif ($consensus eq "falcon") {
      if ($layCounter == 0) {
         open(CNS, "> $prefix.cns.in") or die("Couldn't open '$prefix.cns.in'", undef);
      } else {
         open(CNS, ">> $prefix.cns.in") or die("Couldn't open '$prefix.cns.in'", undef);
      }
      $id =~ s/^\s+|\s+$//g;
      $cns =~ s/^\s+|\s+$//g;
      print CNS "$id $cns\n";
      for (keys %hash) {
         my $str = $hash{$_};
         my @tokenized = split('\s+', $str);
         my $aln =  $tokenized[2];
         $aln =~ s/^\s+|\s+$//g;

         print CNS "$_ $aln\n";
      }
      print CNS "+ +\n";
      close(CNS);
      $layCounter++;
      if ($layCounter > getGlobal("batch")) {
         flushFalcon($output);
         open(CNS, "> $prefix.cns.in") or die("Couldn't open '$prefix.cns.in'", undef);
         close(CNS);
         $layCounter = 0;
      }
   }
}
 
sub processTig($$$) {
   my $input = shift @_;
   my $seq = shift @_;
   my $output = shift @_;

   my %seqs = {};
   my $id = 0;
   my $fastaSeq = "";
   open(F, "< $seq") or die("Couldn't open '$seq'", undef);
   while (<F>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);

      if (m/\>/) {
         if ($fastaSeq ne "") {
            $seqs{$id} = $fastaSeq;
         }
         $fastaSeq = "";
         $id = (split '\s+')[0];
         $id =~ s/>//g;
         if ($id =~ m/,/) {
            $id = (split ',', $id)[1]; 
         } elsif ($id =~ m/utg/) {
            $id =~ s/utg//g;
         }
      } else {
         $fastaSeq = $fastaSeq . $_;
      }      
   }
   if ($fastaSeq ne "") {
      $seqs{$id} = $fastaSeq;
   }
   close(F);
   my $count = (scalar keys %seqs);
   print STDERR "Loaded " . $count . " sequences\n";

   # first loop to get sizes
   my %lens = {};
   my $utgID = 0;
   my $lastEnd = 0;
   my $isUTG = 0;
   open(F, "< $input") or die("Couldn't open '$input'", undef);
   while (<F>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);
      my @tokenized =  split '\s+';
      if ($tokenized[0] eq "unitig") {
         $isUTG = 1;
         if ($lastEnd != 0) {
            $lens{$utgID} = $lastEnd;
         }
         $utgID = $tokenized[1];
      } elsif ($tokenized[0] eq "contig") {
         $isUTG = 0;
         if ($lastEnd != 0) {
            $lens{$utgID} = $lastEnd;
         }
         $utgID = $tokenized[1];
      } elsif ($tokenized[0] eq "FRG") {
         $lastEnd = ($tokenized[13] > $tokenized[14] ? $tokenized[13] : $tokenized[14]);
      }
   }
   if ($lastEnd != 0) {
      $lens{$utgID} = $lastEnd;
   }
   close(F);

   my $counter = 0;
   open(OUT, "> $output") or die ("Couldn't open '$output'", undef);
   open(F, "< $input") or die("Couldn't open '$input'", undef);
   while (<F>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);
      my @tokenized =  split '\s+';
      if ($tokenized[0] eq "unitig") {
         $lastEnd = $lens{$tokenized[1]};
         $utgID = $tokenized[1];
         $isUTG = 1;
         $counter++;
      } elsif ($tokenized[0] eq "contig") {
         $lastEnd = $lens{$tokenized[1]};
         $utgID = $tokenized[1];
         $isUTG = 0;
         $counter++;
      } elsif ($tokenized[0] eq "FRG" && $isUTG != 0) {
           my $min = ($tokenized[13] < $tokenized[14] ? $tokenized[13] : $tokenized[14]);
           my $max = ($tokenized[13] > $tokenized[14] ? $tokenized[13] : $tokenized[14]);
           my $fwd = ($tokenized[13] < $tokenized[14] ? "+" : "-");
           my $seq = $seqs{$tokenized[4]};
           if ($fwd eq "-") {
              $seq =~ tr/ATGCatgcNn/TACGtacgNn/; 
              $seq = reverse $seq;
              $fwd = "+";
           } 
           print OUT "$tokenized[4] utg_$utgID $fwd $lastEnd $min $max $seq $seq\n";
      } elsif ($tokenized[0] eq "UTG" && $isUTG == 0) {
           my $min = ($tokenized[6] < $tokenized[7] ? $tokenized[6] : $tokenized[7]);
           my $max = ($tokenized[6] > $tokenized[7] ? $tokenized[6] : $tokenized[7]);
           my $fwd = ($tokenized[6] < $tokenized[7] ? "+" : "-");
           my $seq = $seqs{$tokenized[4]};
           if ($fwd eq "-") {
              $seq =~ tr/ATGCatgcNn/TACGtacgNn/;
              $seq = reverse $seq;
              $fwd = "+";
           }
           print OUT "$tokenized[4] utg_$utgID $fwd $lastEnd $min $max $seq $seq\n";
      }
   }
   close(F);
   close(OUT);
   print STDERR "Converted $counter unitigs\n";
}

sub processLayouts($$$) {
   my $input = shift @_;
   my $output = shift @_;
   my $sensitive = shift @_;

   my $coverage = getGlobal("maxCoverage");
   my $consensus = getGlobal("consensus");
   my $lastEnd = 0;
   my $lastOffset=0;
   my $cns = "";
   my $lastID = "";
   my %toAlign = {};
   my %bpCov = {};

   my $F;
   if (defined($input)) {
      open($F, "< $input") or die("Couldn't open '$input'", undef);
   } else {
      $F = *STDIN;
   }
   while (<$F>) {
      s/^\s+//;
      s/\s+$//;

      next if (m/^\s*\#/);
      next if (m/^\s*$/);
      my @tokenized =  split '\s+';
      my $hashLen = scalar keys %toAlign;

      if ($lastID eq "" || $lastID ne $tokenized[1]) {
         if ($lastID ne "") { processLayout($output, $lastID, $cns, \%toAlign, $sensitive); } 
         for (keys %toAlign) { delete $toAlign{$_}; }
         $cns = ' ' x ($tokenized[3] + 1);
         for (keys %bpCov) { delete $bpCov{$_}; }
         $lastEnd = 0;
         $lastOffset = 0;
      }
      my $end = $tokenized[5];
      my $start = $tokenized[4];
      if ($start < 0) { 
         $lastOffset = -1*$start;
         $cns = ' ' x ($tokenized[3] + 1 + $lastOffset);
      }
      if ($end >= $lastEnd) {
         if ($start > $lastEnd) {
            $start = $lastEnd;
         }
         if ($end - $start > length($tokenized[7])) {
            $end = $start + length($tokenized[7]);
         }
         $lastEnd = $end;
         substr($cns, $start+$lastOffset, $end-$start) = $tokenized[7];
      }
      my $belowCov = 0;
      if ($sensitive && $consensus eq "pbdagcon") {
         for (my $i = $start; $i < $end; $i++) {
            if (!defined($bpCov{$i})) {
               $bpCov{$i} = 0;
            }
            if ($bpCov{$i} < $coverage) {
               $belowCov++;
            }
            $bpCov{$i} = $bpCov{$i} + 1;
         }
      }
      if ($belowCov > length($tokenized[6])*0.2 || $consensus ne "pbdagcon" || $sensitive == 0) {
         $toAlign{$tokenized[0]} = $start+$lastOffset . " " . $end . " " .  $tokenized[6];
      }
      $lastID = $tokenized[1];
   }
   if ($lastID ne "") {
      processLayout($output, $lastID, $cns, \%toAlign, $sensitive);
   }
   flushFalcon($output);
   flushPBDAGCON($output);
   flushPBUTGCNS($output);
   if (defined($input)) {
      close($F);
   }
}

my $err = 0;
setGlobal("consensus", "pbdagcon");
setGlobal("prefix", "tmp");
setGlobal("coverage", 0);
setGlobal("threads", 1);
setGlobal("maxCoverage", 5000);

my $input = undef;
my $output = undef;
my $inputSequences = undef;

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

    } elsif ($arg eq "-pythonPath") {
        setGlobal("falconpath", shift @ARGV);

    } elsif ($arg eq "-input") {
       $input = shift @ARGV;

    } elsif ($arg eq "-output") {
        $output = shift @ARGV;

    } elsif ($arg eq "-sequence") {
        $inputSequences = shift @ARGV;

    } elsif ($arg eq "-consensus") {
        setGlobal("consensus", shift @ARGV);

    } else {
       print STDERR "Unknown parameter " + $err + "\n";
       $err++;
    }
 }

if ($err > 0 || !defined($output)) {
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

my $consensus = getGlobal("consensus");
if ($consensus eq "falcon") {
   my $falconPath = getGlobal("falconpath");
   if (! -e "$falconPath/falcon_sense.py") {
       print STDERR "Error: selected falcon but it was not available, using pbdagcon instead\n";
       setGlobal("consensus", "pbdagcon");
   }
}
if ($consensus eq "pbutgcns") {
   my $path = getGlobal("smrtpath"); 
   if (!defined($path) || $path eq "" || ! -e "$path/pbutgcns") {
      # try to use path
      my $binPath = `which pbutgcns`;
      chomp $binPath;
      my @t = split '/', "$binPath";
      pop @t;                      #  pbutgncs
      $path = join '/', @t;  #  path to the assembler
   }
   if (! -e "$path/pbutgcns") {
       print STDERR "Error: selected pbutgcns but it was not available, using pbdagcon instead\n";
       setGlobal("consensus", "pbdagcon");
   } else {
      setGlobal("smrtpath", $path);
   }
}
if ($consensus eq "pbdagcon") {
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
   if ((! -e "$BLASR/blasr" || ! -e "$BLASR/pbdagcon")) {
      print STDERR "Error, could not find blasr/smrtportal in path. Please update your path or specify a location using -path\n";
      exit(1);
   }
}

# clear output file
open(O, ">", $output);
close(O);

# when we were given sequences in a fasta file, assume we are processing a tig store layout, otherwise assume we have a layout already
if (defined($inputSequences)) {
   my $out = getGlobal("prefix") . ".layout";
   processTig($input, $inputSequences, $out);
   processLayouts($out, $output, 0);
} else {
   processLayouts($input, $output, 1);
}
# cleanup
system("rm -f " . getGlobal("prefix") . "*");
