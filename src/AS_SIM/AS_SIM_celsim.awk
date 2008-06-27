BEGIN {
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
    notdone = 1;
    srand();
    keepDNA = 0; # Don't nuke the dna file
    fileComment = "(none)"
    sawfileseed = 0;# Have we seen an (optional) .seed section?
    sawcmdlineseed = 0; # Have we seen a command line seed
    seed    = int(1000*rand());
    notations = 0;
    comment = 1;
    fastaOutput = 0;
    uniform = 0;
    base    = 1;
    prompt = 0;
    protoflag = " -P ";
    split("$Revision: 1.5 $", revFields, " ");
    revision = revFields[2];
    split("$Date: 2008-06-27 06:29:21 $", dateFields, " ");
    date = dateFields[2];
    print STDERR "date = " date " revision = " revision;
    while (base < ARGC-1)
      { if (ARGV[base] ~ "-F")
	comment = 0;
      else if (ARGV[base] ~ "-s" && base+1 < ARGC-1)
	{ cmdlineseed = ARGV[base+1];
	sawcmdlineseed = 1;
	base += 1;
	}
      else if (ARGV[base] ~ "-d"){
	print "* Saving DNA File\n";
	keepDNA = 1;
      }
      else if (ARGV[base] ~ "-p"){
	print "* Will Prompt before clobbering files\n";
	prompt = 1;
      }
      else if (ARGV[base] ~ "-u"){
	print "* Uniform Distance Distribution (Normal is default)\n";
	uniform = 1;
      }
      else if (ARGV[base] ~ "-f"){
	print "* FASTA output (Assembler ProtoI/O is default)\n";
	fastaOutput = 1;
	protoflag = " ";
      }
      else if (ARGV[base] ~ "-n"){
	print "* Normal Distance Distribution (default)\n";
	uniform = 0;
      }
      else if (ARGV[base] ~ "-b" && base+1 < ARGC-1)
	{ cmdlinebatchsize = ARGV[base+1];
	sawcmdlinebatchflag = 1;
	base += 1;
	}
      else
      { print " *** Illegal argument '%s'", ARGV[base];
      notdone = 0;
      exit (1);
      }
        base += 1;
      }
    if (base != ARGC-1)
      { print " *** Usage: celsim [-F] [-s #] [-d] [-p] [-u] [-f] [-n] [-b <batch size>] file";
        notdone = 0;
        exit (1);
      }

    pfile = ARGV[1] = ARGV[ARGC-1];
    ARGC  = 2;
    numFields = split(pfile,fields,".");
    if(numFields == 1){   # No Extension
      tryfile = pfile ".sim";
      prefix = pfile;
    }
    else{
      # print "file extension is " fields[numFields];
#      if("spec" != fields[numFields] &&
	# "sim" != fields[numFields]){
n#	print "File extension " fields[numFields] " is invalid -- try .sim\n";
#	exit(1);
#      }else{

	prefix = "";
	for(i = 1; i < numFields; i++){
	  prefix = prefix fields[i];
	}

      }
#    }


    frgfile = prefix ".frg";
    outfile = prefix ".cms";
#    notfile = prefix ".not." seed;
    notfile = prefix ".not";
    adlfile = prefix ".adl";
#    comfile = prefix ".cms." seed;
    comfile = prefix ".com";
#    qltfile = prefix ".qlt." seed;
    qltfile = prefix ".qlt"
#CMM    dnafile = ".dna." seed;
    dnafile = prefix ".dna"

    phase = -1;
    npoly = 0;
    ptotl = 0.;
    nlibs = 0;
    acnum = 1;

  }

/^\.seed/ && notdone {
    if (phase != -1)
      { print " *** Seed specification must come first.";
        cleanup(1);
      }
    # Do not change phase...we are still before .dna
    # print "SawSeed!\n";
    sawfileseed = 1;
  }
/^\.comment/ && notdone {
    fileComment = "";
    for (i = 2; i <= NF; i++)
      fileComment = fileComment " " $i;
    print "FileComment" fileComment;
  }

/^\.dna/ && notdone {
    if (phase != -1)
      { print " *** Dna specification must come first.";
        cleanup(1);
      }
    phase = 0;
    # Cmd line uber alles
    if(sawcmdlineseed == 1){
      seed = cmdlineseed;
      print "Command line seed won: seed = " seed;
    } else if (sawfileseed == 2){ #file .seed section over rand
      seed = fileseed;
      # print "File seed won: seed = " seed;
    }
}

/^\.poly([\t ]+[0-9.]+)+$/ && notdone {
    if (phase == -1)
      { print " *** Missing dna specification.";
        cleanup(1);
      }
    if (phase == 2)
      { print " *** Poly specification after sample spec.";
        cleanup(1);
      }
    if (phase == 0)
      make_dna();
    else
      make_polys(pcopy);
    phase = 1;
    for (i = 2; i <= NF; i++)
      pweight[npoly+(i-1)] = $i;
    pcopy = NF-1;
    lines = 0;
  }

/^\.sample$/ && notdone {
    if (phase == -1)
      { print " *** Missing dna and polymorphism specifications.";
        cleanup(1);
      }
    if (phase == 0)
      { make_dna();
        npoly = 1;
        ptotl = 1.;
        pweight[1] = 1.;
	if(notations > 0)
	  close(notfile);
        print "# NO POLYMORPHISMS APPLIED, USED BASE SEQUENCE" >>outfile;
        print "#" >>outfile;
	cmd = "cp " dnafile  " " prefix ".poly.1";
	print "Executing: " cmd;
        system(cmd);
      }
    else if (phase == 1)
      make_polys(pcopy);
    else{
      print "* Invoking make_sample from .sample";
      make_sample();
    }
    phase = 2;
    items = 0;
    bacends = 0;
  }

/^\.bacends$/ && notdone {
    if (phase == -1)
      { print " *** Missing dna and polymorphism specifications.";
        cleanup(1);
      }
    if (phase == 0)
      { make_dna();
        npoly = 1;
        ptotl = 1.;
        pweight[1] = 1.;
	if(notations > 0)
	  close(notfile);
        print "# NO POLYMORPHISMS APPLIED, USED BASE SEQUENCE" >>outfile;
        print "#" >>outfile;
	cmd = "cp " dnafile " "  prefix ".poly.1";
	print "Executing: " cmd;
        system(cmd);
      }
    else if (phase == 1)
      make_polys(pcopy);
    else{
      # print "* Invoking make_sample from .sample";
      make_sample();
    }
    phase = 2;
    items = 0;
    bacends = 1;
  }

END {
  if (notdone) {
    if (phase != 2) {
      print " *** Generated no data! (phase = " $phase ")";
      cleanup(1);
    }
    make_sample();
    if(npoly == 1 && notations > 0 && !fastaOutput){
      close(sname);  # close the file, so we can annotate it
      perlcmd = "perl \$AS_ROOT/bin/AS_SIM_labelfrag.pl " outfile " " sname " " notfile;
      if(system(perlcmd) != 0){
	print "*** Error in fragment labeling\n" >> outfile;
	cleanup(1);
	exit(1);
      }
    }else if (notations > 0){

      if(!fastaOutput)
	print "*** Poly > 1 -- No fragment labelling applied\n";
    }

    if(!fastaOutput && sawcmdlinebatchflag){
      system("tobatches -b " cmdlinebatchsize " " prefix);
    }
    cleanup(0);
  }
}

/^([^.]|\.[0-9])/ && notdone {
# Process a section body .......
  if (phase == -1 && sawfileseed == 1){
    fileseed = $1;
    sawfileseed = 2; # Read only one seed in this section
    # print "Seed = " seed "\n";
  }else  if (phase == 0){
    if($1 == "@"){
      notations++; # track number of notations
      print $2 > notfile;
      $1 = "";  # Strip the @ sign, so frag doesn't see it
    }
      print > comfile
  } else if (phase == 1)
      { lines += 1;
        if (phase == 1 && lines == 1) print dnafile >comfile;
        print >comfile;
      }
    else if (phase == 2)
      { if (items == 0 && NF >= 1)
          {  j = 2;
	  fnum = $1;
	  items++;
#	  print STDERR "items " items " fnum = " fnum ;
	    }
        else
        j = 1;
        for (i = j; i <= NF; i++){
#	  print STDERR "i = " i " NF = " NF " items = " items " text = " $i;
	  if(items < 4 || items > 7 ){
	    fragspec = fragspec " " $i;
	  }else{
	    massageSpec = massageSpec " " $i;
	  }
	  items++;
	}
      }
  }

function make_dna()
{ close(comfile);
  print " +++ Building DNA sequence";
  if (system("$AS_ROOT/bin/frag -s " seed protoflag " " comfile " >" dnafile) != 0)
    { print " *** Dna spec error";
      cleanup(1);
    }
    if(uniform)
      uni = " (uniform) ";
    else
      uni = " (gaussian) ";
    if(!fastaOutput){
      print "Celsim " revision uni " " date > adlfile;
      perlcmd = "perl \$AS_ROOT/bin/AS_SIM_extractLength.pl  < " dnafile " >> " adlfile ;
      print fileComment >> adlfile;
# print STDERR perlcmd;
      if(system(perlcmd) != 0){
	print " *** Error in extractLength.pl";
	cleanup(1);
      }
     # print "*** Celsim Input was:\n" >> adlfile;
      perlcmd = "echo \" ********Celsim Input ******** \">> "adlfile;
     #      print STDERR perlcmd;
      if(system(perlcmd) != 0){
	print STDERR " *** Error outputing comment";
	cleanup(1);
      }
      perlcmd = "cat " pfile " >>" adlfile ;
      if(system(perlcmd) != 0){
	print STDERR " *** Error in cat input to adlfile";
	cleanup(1);
      }
      perlcmd = "echo \" ********End Celsim Input ******** \">> "adlfile;
      if(system(perlcmd) != 0){
	print STDERR " *** Error outputting comment";
	cleanup(1);
      }


    }
  if (comment)
    { if (overwrite(outfile)) return;
      print "# Celsim V" revision " " fileComment " " >> outfile
      print "# DNA SEQUENCE GENERATION" >>outfile;
      system("awk '/^#/ && ($0 !~ /^# DNA Sequence:/)' " dnafile " >>" outfile);
    }
}

function make_polys(pcopy)
{ close(comfile);
  for (i = 1; i <= pcopy; i++)
    { npoly += 1;
      ptotl += pweight[npoly];
      print " +++ Making polymorphism " npoly;
      command = "$AS_ROOT/bin/poly -s " seed+npoly " " comfile " > " prefix".poly." npoly ;
      v = system(command);
      if (v != 0)
        { print STDERR " Error Invoking " command " ...exiting";
          cleanup(1);
        }
      if (comment)
        { if (npoly != 1) print "#" >>outfile;
          print "# POLYMORPHISM GENERATION: WEIGHT = " pweight[npoly] >>outfile;
          system("awk '/^#/' " prefx ".poly." npoly " >>" outfile);
        }
    }
}

# massageSpec
#   1  min prefix unclear range
#   2  max prefix unclear range
#   3  min suffix unclear range
#   4  max suffix unclear range

# fragspec gets split as follows:
#   1  min frag len
#   2  max frag len
#   3  F/R odds
#   4  Error ramp Min
#   5  Error ramp Max
#   6  Ins odds
#   7  Del odds
#   8  Single Odds
#   9  min insert length
#   10 max insert length
#   11 False mate rate
#
function make_sample()
{ nlibs += 1;
  sname = frgfile ;
  if (nlibs == 1)
    if (overwrite(sname))
      return;

  if(nlibs == 1 && !fastaOutput){
    command = "$AS_ROOT/bin/outputADT < " adlfile "  >" sname;
    print STDERR "adlfile = " adlfile " sname= " sname " command = " command;
      v = system(command);
      if (v != 0)
        { print STDERR " Error Invoking " command " ...exiting";
          cleanup(1);
        }
  }
  {
    split(fragspec,splitFragSpec);
    split(massageSpec,splitMassageSpec);

    med = int((splitFragSpec[10] + splitFragSpec[9])/2);
    dta = int((splitFragSpec[10] - splitFragSpec[9])/2);

    if(items > 12 && !fastaOutput){
      command = "$AS_ROOT/bin/outputDST " acnum " " med " " dta  " >>"  sname;

# print "fragspec " fragspec " \nend of fragpsec\n";
      v = system(command);
      if (v != 0)
        { print STDERR " Error Invoking " command " ...exiting";
	cleanup(1);
        }
      distid = acnum;
      acnum += 1;
    }
  }

  cwght = 0.;
  lind  = 0;
  for (i = 1; i <= npoly; i++)
    {
      close(comfile);
      print "< " prefix ".poly." i  >comfile;
      print STDERR "< " prefix ".poly." i ;
      cwght += pweight[i];
      nind   = int(fnum*cwght/ptotl);
      nfrag[i] = nind-lind;
      print nfrag[i] >comfile;
      print STDERR nfrag[i] ;
      lind   = nind;
      print fragspec >comfile;
      print STDERR fragspec;
      close(comfile);
      c = seed + nlibs*npoly + i;
# frag could have a command line argument with a label
# for the particular polymorphism we're generating from.
# This would facilitate labeling the fragments with their poly
# origin.
      if(uniform)
	uni = " -u ";
      else
	uni = " -n ";
      if(bacends)
	bac = " -b ";
      else
	bac = " ";
      # print "### uni = " uni;
      if(!fastaOutput){
	command = "$AS_ROOT/bin/frag -s " c protoFlag uni " -F -N " comfile " | $AS_ROOT/bin/massage -q " qltfile " ";
	command = command  uni bac distid " " acnum " " massageSpec " >>" sname;
       print STDERR command;
	print STDERR massageSpec;
	print STDERR fragspec;

      }else{
	command = "$AS_ROOT/bin/frag -s " c uni " -F -N " comfile " >> " sname
      }
     print " +++ Generating library " nlibs " from poly " i;
      v = system(command);
      acnum += nfrag[i];
      if (v != 0)
        { print STDERR " Error Invoking " command " ...exiting";
          cleanup(1);
        }
    }
  if (comment)
    { printf "#\n" >>outfile;
      printf "# FRAGMENT LIBRARY %d\n", nlibs >>outfile;
      printf "#\n" >>outfile;
      printf "# Fragments   Seed   Poly Source\n" >>outfile;
      printf "# ---------   ----   -----------\n" >>outfile;
      for (i = 1; i <= npoly; i++)
        printf "#  %8d  %5d   %s.poly.%d.%d\n",
               nfrag[i], seed + nlibs*npoly + i, prefix, i, seed  >>outfile;
      printf "# ---------\n#  %8d\n", fnum >>outfile;
      printf "#\n" >>outfile;
      printf "#    Length Range = [%d,%d], F/R odds = %.2f/%.2f\n",
             splitFragSpec[1], splitFragSpec[2], splitFragSpec[3], 1.-splitFragSpec[3] >>outfile;
      printf "#\n" >>outfile;
      printf "# Clear Range Characteristics:\n" >> outfile;
      printf "#    Prefix Unclear [%d,%d]\n", splitMassageSpec[1], splitMassageSpec[2] >> outfile
      printf "#    Suffix Unclear [%d,%d]\n", splitMassageSpec[3], splitMassageSpec[4] >> outfile
      printf "# Edit Characteristics:\n" >>outfile;
      printf "#    Error Ramp = %.2f->%.2f, ", splitFragSpec[4], splitFragSpec[5] >>outfile;
      printf "Ins/Del/Sub Odds = %.2f/%.2f/%.2f\n",
             splitFragSpec[6], splitFragSpec[7], 1.-(splitFragSpec[6]+splitFragSpec[7]) >>outfile
      printf "#\n" >>outfile;
      if (items > 12)
        { printf "# Dual-End Inserts:\n" >>outfile;
          printf "#    Single Odds = %.2f\n", splitFragSpec[8] >>outfile;
          printf "#    Insert Range = [%d,%d]\n", splitFragSpec[9], splitFragSpec[10] >>outfile;
          printf "#    Pairing Error Rate = %.2f\n", splitFragSpec[11] >>outfile;
          printf "#\n" >>outfile;
        }
    }
  fragspec = "";
  massageSpec = "";
}

function cleanup(ecode)
{
  exit(ecode); # TEMPORARY, for debug
  if (phase >= 0){
    if(keepDNA == 0)
      system("rm -f " dnafile);
    if(notations > 0)
      system("rm " notfile);
  }
  if (phase >= 1)
    for (i = 1; i <= npoly; i++)
      system("rm " prefix ".poly." i );
  system("rm " comfile);
  system("rm " qltfile);
#  if(!fastaOutput)
#    system("rm " adlfile);
  notdone = 0;
  exit (ecode);
}

function overwrite(name)
{ if (system("test -e " name) == 0)
    {
      if(prompt){
	printf "  Overwrite %s [y/n]? ", name;
	getline answer <"/dev/stdin";
	if (answer ~ "^y$")
	  system("rm " name);
	else
	  return (1);
      }else{
	  system("rm " name);
      }
    }
  return (0);
}
