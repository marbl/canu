#!/usr/bin/perl 
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
#********************************************************************
#        Script:  
#   Description:  Perl script that generates and runs a simfile
#                 for different parameter and outputs the results
#                 in a gnuplot file
#    Programmer:  Knut Reinert
#       Written:  10 June 99
#  Last Revised:  10 June 99
#********************************************************************

$file = shift(@ARGV);
$screen = shift(@ARGV);


$num_repeats = 25;
$repeat_size = 2000;
$random_sequence_size = 110;
$genome_size = 500000;

for($coverage=5; $coverage<=20; $coverage+=5){
  open(GNUDATAHANDLE,"> $file.gpd.$coverage");
  open(GNU3COMMANDHANDLE,"> $file.gpc.$coverage");
  if( $screen ne "s" ){
    print GNU3COMMANDHANDLE "set output \"$file.$coverage.ps\"\n";
    print GNU3COMMANDHANDLE "set terminal postscript \n";
  }
  print GNU3COMMANDHANDLE "set xlabel \"mutation error\"\n";
  print GNU3COMMANDHANDLE "set ylabel \"sequencing error\"\n";
  print GNU3COMMANDHANDLE "set zlabel \"number of chunks\"\n";
  print GNU3COMMANDHANDLE "set label \"length=$genome_size\" at screen 0.02,0.97\n";
  print GNU3COMMANDHANDLE "set label \"num_rep=$num_repeats\" at screen 0.02,0.94 \n";
  print GNU3COMMANDHANDLE "set label \"rep_length=$repeat_size\" at screen 0.02,0.91\n";
  print GNU3COMMANDHANDLE "set label \"rand_seq_length=$random_sequence_size\" at screen 0.02,0.88\n";
  print GNU3COMMANDHANDLE "splot \"$file.gpd.$coverage\" w lines\n";
  if( $screen eq "s" ){
    print GNU3COMMANDHANDLE "pause -1;\n";
  }
  close GNU3COMMANDHANDLE;

  for($seq_err=0.005; $seq_err <= 0.02; $seq_err += 0.005){
    open(TWODHANDLE,"> $file.gpd.$coverage.$seq_err");

    open(GNU2COMMANDHANDLE,"> $file.gpc.$coverage.$seq_err");
    if( $screen ne "s" ){
      print GNU2COMMANDHANDLE "set output \"$file.$coverage.$seq_err.ps\"\n";
      print GNU2COMMANDHANDLE "set terminal postscript \n";
    }
    print GNU2COMMANDHANDLE "set xlabel \"mutation error\"\n";
    print GNU2COMMANDHANDLE "set ylabel \"number of chunks\"\n";
    print GNU2COMMANDHANDLE "set label \"length=$genome_size\" at screen 0.75,0.85\n";
    print GNU2COMMANDHANDLE "set label \"num_rep=$num_repeats\" at screen 0.75,0.82 \n";
    print GNU2COMMANDHANDLE "set label \"rep_length=$repeat_size\" at screen 0.75,0.79\n";
    print GNU2COMMANDHANDLE "set label \"rand_seq_length=$random_sequence_size\" at screen 0.75,0.76\n";
    print GNU2COMMANDHANDLE "plot \"$file.gpd.$coverage.$seq_err\" w lines\n";
    if( $screen eq "s" ){
      print GNU2COMMANDHANDLE "pause -1;\n";
    }
    close GNU2COMMANDHANDLE;

    for($mut_err=0.0; $mut_err <= 0.2; $mut_err += 0.005){

      $out = `create-rep-sim $file $coverage $seq_err $mut_err $num_repeats $repeat_size $random_sequence_size $genome_size`;

      print $out,"\n";

      open(INPUT,"$file.cga");

      while( <INPUT> )
	{
	  chop;
	  if( $_ =~ /Total number of chunks/) {
	    ($nc,$rest) = split;
	  }
	}
      close INPUT;
      print "$mut_err $seq_err $nc \n";
      flock (TWODHANDLE, 2);
      print TWODHANDLE "$mut_err $nc \n";
      flock (TWODHANDLE, 8);
      flock (GNUDATAHANDLE, 2);
      print GNUDATAHANDLE "$mut_err $seq_err $nc \n";
      flock (GNUDATAHANDLE, 8);

      # delete all files (in case a command does not succed)
#      `rm -rf $file-GKP`;
#      `rm -rf $file-OVL`;
#      `rm  $file.cam`;
#      `rm  $file.cga`;
#      `rm  $file.cgb`;
#      `rm  $file.cms`;
#      `rm  $file.fgb`;
#      `rm  $file.frg`;
#      `rm  $file.inp`;
#      `rm  $file.ovl`;
#    `rm  $file.sim`;
    }
    close TWODHANDLE;
    flock (GNUDATAHANDLE, 2);
    print GNUDATAHANDLE "\n";
    flock (GNUDATAHANDLE, 8);
  }
  close GNUDATAHANDLE;

}
    

