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

# scan commandline
use Getopt::Std;
getopts('c:f:') or die "Wrong argument list : call split_celamy.pl -c \# -f input_celamy\n";
$opt_c or  die "Wrong argument list : call split_celamy.pl -c \# -f input_celamy\n";
$opt_f or  die "Wrong argument list : call split_celamy.pl -c \# -f input_celamy\n";
use Cwd;
$cwd = cwd();
$contigs_per   = $opt_c;
$input_file   = $opt_f;

$basename = $input_file;
$basename =~ s/[.].*//;
print "Basename of file is: $basename\n";
# first, split input into header, body and links
unlink("head.t");
unlink("body.t");
unlink("links.t");
$split_cmd = "gawk '{if (match(\$2,/^C/) ) { print >> \"head.t\";next;} if (match(\$1,/LNK/) ){print >> \"links.t\";next;} print >> \"body.t\" }' $opt_f";
#print "$split_cmd\n";
system($split_cmd);


$contig_line_cmd =  "gawk 'BEGIN { contigs=0 } { contigs += NF-1;if (contigs > $contigs_per) { print NR; contigs = 0 } }' links.t";
#print "$contig_line_cmd\n";
open CLINES, "$contig_line_cmd |";

@clines = <CLINES>;

#print "Arg_count on clines is $#clines\n";
print @clines;

system("rm links[0-9]*.t");
system("rm body[0-9]*.t");

$item=0;
$cline_cmd = "";
foreach $cline (@clines) {
  $cline =~ s/\n//;
  $cline_cmd .= "if (NR < $cline ) { print >> \"links${item}.t\";next;}";
  $item++;
}
$cline_cmd = "gawk '{ $cline_cmd print >> \"links${item}.t\" }' links.t";

system($cline_cmd);

$filecount = $item+1;
print "filecount is $filecount\n";

@first_ScafCtg;
@coord;
@lineno;
$item =0;
$coord_cmd = "";
for ($i=0;$i<$filecount;$i++) {
   open LINKFILE, "head -1 links$i.t |";
   $link = <LINKFILE>;
   if ( $link =~ /LNK:/ ) { 
     @pline = split ' ',$link;
     shift(@pline);
     $first_ScafCtg[$item] = shift(@pline); 
     $coord_cmd .= "if ( match(\$1,/^$first_ScafCtg[$item]:/) ) { print NR,\":\",\$0; next;}";
     $item++;
   }
}
$coord_cmd = "gawk '{ $coord_cmd }' body.t";
#print "Coord command: ",$coord_cmd;
open COORDS, "$coord_cmd |";
$item = 0;
while ( <COORDS> ) {
     $cline = $_;
     #print "Contig line:",$cline;
     @pcline = split ':',$cline;
     $lineno[$item] = shift @pcline;
     shift @pcline;
     @ppcline = split ' ',shift(@pcline);
     $coord[$item] = shift @ppcline;
     #print "Lineno: ",$lineno[$item]," Coord: ",$coord[$item],"\n";
     $item++;
}

#print "item is ",$item,"\n";

for ( $lnk=1;$lnk<$item;$lnk++ ) {
  $bodyfile = sprintf("body%d.t",$lnk -1);
  $c = $coord[$lnk-1];
  $cmd .= "if ( NR<$lineno[$lnk] ) {\$2-=$c;\$4-=$c;print >> \"$bodyfile\";next;}"
}

$bodyfile = sprintf("body%d.t",$item-1);
$c = $coord[$item-1];
$cmd .= " \$2-= $c ;\$4-= $c ;print >> \"$bodyfile\"; ";

$cmd = "gawk '{ $cmd }' body.t";
print "$cmd\n";
system($cmd);

for ( $lnk=0;$lnk<$item;$lnk++) {
  system("cat head.t body$lnk.t links$lnk.t >> ${basename}_$lnk.cam");
  unlink("body$lnk.t");
  unlink("links$lnk.t");
}

