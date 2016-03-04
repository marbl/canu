
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_MER/mercy.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2007-MAR-19 to 2014-APR-11
 *      are Copyright 2007-2009,2011,2013-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AS_global.H"

#include "libmeryl.H"

//  The categories depend on the type of input (fragments or contigs):
//
//  0 -- no count, mer not present
//  1 -- single copy
//  2 --   2 ->  10 copies (contigs) --    2 ->   2mode copies (frags)
//  3 --  11 -> 100 copies (contigs) --      ->  10mode copies (frags)
//  4 -- 101+ copies (contigs)       --      -> 100mode copies (frags)
//  5                                --      -> infinity copies (frags)

//  You'll also need to modify compare() and output() if you change this.
#define NUMCATEGORIES  6

//  The output files are global for convenience.  Otherwise, we'd be passing
//  them to compare() for every single mer.
//
bool   dumpFlag = false;
FILE  *dumpSCZF = 0L;
FILE  *dumpMCZF = 0L;
FILE  *dumpMCSF = 0L;
FILE  *dumpMCMF = 0L;
char   merstring[1024];

uint32
findMode(char *name) {
  merylStreamReader  *M = new merylStreamReader(name);
  uint32             *H = new uint32 [16384];

  fprintf(stderr, "Finding mode of '%s'\n", name);

  for (uint32 i=0; i<16384; i++)
    H[i] = 0;

  while (M->validMer()) {
    if (M->theCount() < 16384)
      H[M->theCount()]++;
    M->nextMer();
  }

  uint32  mi = 16;
  for (uint32 i=mi; i<16384; i++)
    if (H[i] > H[mi])
      mi = i;

  fprintf(stderr, "Mode of '%s' is "F_U32"\n", name, mi);

  return(mi);
}


void
compare(merylStreamReader *F,
        merylStreamReader *C,
        kMer              &minmer,
        uint32             mode,
        uint32             R[NUMCATEGORIES][NUMCATEGORIES]) {
  uint32  Ftype = 0;
  uint32  Ctype = 0;
  kMer    Fmer  = F->theFMer();
  kMer    Cmer  = C->theFMer();
  uint32  Fcnt  = F->theCount();
  uint32  Ccnt  = C->theCount();

  if (Fcnt == 0)
    Ftype = 0;
  else if (Fcnt == 1)
    Ftype = 1;
  else if (Fcnt <= 2*mode)
    Ftype = 2;
  else if (Fcnt <= 10*mode)
    Ftype = 3;
  else if (Fcnt <= 100*mode)
    Ftype = 4;
  else
    Ftype = 5;

  if (Ccnt == 0)
    Ctype = 0;
  else if (Ccnt == 1)
    Ctype = 1;
  else if (Ccnt <= 10)
    Ctype = 2;
  else if (Ccnt <= 100)
    Ctype = 3;
  else
    Ctype = 4;

  //  If the mer isn't valid, we hit the end of the file, and the mer
  //  thus (obviously) isn't in the file.
  //
  if (F->validMer() == false)
    Ftype = 0;
  if (C->validMer() == false)
    Ctype = 0;

  //  If either type is 0, we're done, but only increment the count if
  //  this mer is the minmer.
  //
  if ((Ftype == 0) || (Ctype == 0)) {
    if (((Ftype == 0) && (Cmer == minmer)) ||
        ((Ctype == 0) && (Fmer == minmer))) {
      R[Ftype][Ctype]++;

      //  Save the mer if it's in contigs, but not fragments.
      if (dumpFlag)
        if (Ftype == 0)
          if (Ctype == 1)
            fprintf(dumpSCZF, ">"F_U32"\n%s\n", Ccnt, Cmer.merToString(merstring));
          else
            fprintf(dumpMCZF, ">"F_U32"\n%s\n", Ccnt, Cmer.merToString(merstring));
    }
    return;
  }

  //  If the mers don't agree, we're also done.  If either is the
  //  minmer, note that we saw it.
  //
  if (Fmer != Cmer) {
    if (Fmer == minmer)
      R[Ftype][0]++;
    if (Cmer == minmer) {
      R[0][Ctype]++;

      //  Again, save the mer since it's in contigs, but not fragments.
      if (dumpFlag)
        if (Ctype == 1)
          fprintf(dumpSCZF, ">"F_U32"\n%s\n", Ccnt, Cmer.merToString(merstring));
        else
          fprintf(dumpMCZF, ">"F_U32"\n%s\n", Ccnt, Cmer.merToString(merstring));
    }

    return;
  }

  //  If we're not the minmer, we're done.
  if (Fmer != minmer)
    return;

  //  Otherwise, the mers are in both inputs
  R[Ftype][Ctype]++;

  //  Save the mer if it's in contigs "more" than if in fragments.
  if (dumpFlag) {
    if (Ftype < Ctype)
      if (Ctype == 2)
        fprintf(dumpMCSF, ">"F_U32"\n%s\n", Ccnt, Cmer.merToString(merstring));
      else
        fprintf(dumpMCMF, ">"F_U32"\n%s\n", Ccnt, Cmer.merToString(merstring));

    if ((Ftype == 0) && (Ctype == 1))
      fprintf(dumpSCZF, ">"F_U32"\n%s\n", Ccnt, Cmer.merToString(merstring));
  }
}


void
output(char              *title,
       uint32             mode,
       uint32             R[NUMCATEGORIES][NUMCATEGORIES]) {

  fprintf(stdout, "\n\n%s\n", title);
  fprintf(stdout, "(frags)    |      zero |       one |     <= 10 |    <= 100 |    <= inf | (contigs)\n");

  for (uint32 i=0; i<6; i++) {
    switch (i) {
      case 0:  fprintf(stdout, "zero       ");  break;
      case 1:  fprintf(stdout, "one        ");  break;
      case 2:  fprintf(stdout, "<= 2mode   ");  break;
      case 3:  fprintf(stdout, "<= 10mode  ");  break;
      case 4:  fprintf(stdout, "<= 100mode ");  break;
      case 5:  fprintf(stdout, "<= inf     ");  break;
      default: fprintf(stdout, "?????????  ");  break;
    }
    for (uint32 j=0; j<5; j++)
      fprintf(stdout, "%12"F_U32P, R[i][j]);
    fprintf(stdout, "\n");
  }
}




int
main(int argc, char **argv) {
  merylStreamReader  *AF = 0L;
  merylStreamReader  *TF = 0L;
  merylStreamReader  *AC = 0L;
  merylStreamReader  *DC = 0L;
  merylStreamReader  *CO = 0L;

  uint32              AFmode = 0;
  uint32              TFmode = 0;

  char                dumpSCZFname[1024] = {0};  //  single contig, zero frags
  char                dumpMCZFname[1024] = {0};  //  low contig, zero frags
  char                dumpMCSFname[1024] = {0};  //  medium contig, low frags
  char                dumpMCMFname[1024] = {0};  //  everything else, contig > frags

  bool                beVerbose = false;

  argc = AS_configure(argc, argv);

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-af") == 0) {  //  All frags
      ++arg;
      AFmode = findMode(argv[arg]);
      AF = new merylStreamReader(argv[arg]);
      AF->nextMer();
    } else if (strcmp(argv[arg], "-tf") == 0) {  //  Trimmed frags
      ++arg;
      TFmode = findMode(argv[arg]);
      TF = new merylStreamReader(argv[arg]);
      TF->nextMer();
    } else if (strcmp(argv[arg], "-ac") == 0) {  //  All contigs
      AC = new merylStreamReader(argv[++arg]);
      AC->nextMer();
    } else if (strcmp(argv[arg], "-dc") == 0) {  //  Degenerate contigs
      DC = new merylStreamReader(argv[++arg]);
      DC->nextMer();
    } else if (strcmp(argv[arg], "-co") == 0) {  //  Contigs
      CO = new merylStreamReader(argv[++arg]);
      CO->nextMer();
    } else if (strcmp(argv[arg], "-dump") == 0) {
      arg++;
      dumpFlag = true;
      sprintf(dumpSCZFname, "%s.0.singlecontig.zerofrag.fasta",       argv[arg]);
      sprintf(dumpMCZFname, "%s.1.multiplecontig.zerofrag.fasta",     argv[arg]);
      sprintf(dumpMCSFname, "%s.2.multiplecontig.lowfrag.fasta",      argv[arg]);
      sprintf(dumpMCMFname, "%s.3.multiplecontig.multiplefrag.fasta", argv[arg]);
    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((AF == 0L) && (TF == 0L) && (AC == 0L) && (DC == 0L) && (CO == 0L)) {
    fprintf(stderr, "usage: %s [opts] [-v] [-dump prefix]\n", argv[0]);
    fprintf(stderr, "At least one fragcounts and one contigcounts are needed.\n");
    fprintf(stderr, "          -af | -tf        fragcounts\n");
    fprintf(stderr, "          -ac | -dc | -co  contigcounts \n");
    fprintf(stderr, "Dumping is probably only useful with exactly one frag and\n");
    fprintf(stderr, "one contig, but I'll let you do it with any number.\n");
    exit(1);
  }
  if ((AF == 0L) && (TF == 0L)) {
    fprintf(stderr, "ERROR - need at least one of -af, -tf\n");
    exit(1);
  }
  if ((AC == 0L) && (DC == 0L) && (CO == 0L)) {
    fprintf(stderr, "ERROR - need at least one of -ac, -dc, -co\n");
    exit(1);
  }

  //  Check mersizes.
  //
  uint32  merSize = 0;
  uint32  ms[5] = { 0 };

  if (AF)  merSize = ms[0] = AF->merSize();
  if (TF)  merSize = ms[1] = TF->merSize();
  if (AC)  merSize = ms[2] = AC->merSize();
  if (DC)  merSize = ms[3] = DC->merSize();
  if (CO)  merSize = ms[4] = CO->merSize();

  bool  differ = false;

  if ((ms[0] > 0) && (ms[0] != merSize))  differ = true;
  if ((ms[1] > 0) && (ms[1] != merSize))  differ = true;
  if ((ms[2] > 0) && (ms[2] != merSize))  differ = true;
  if ((ms[3] > 0) && (ms[3] != merSize))  differ = true;
  if ((ms[4] > 0) && (ms[4] != merSize))  differ = true;

  if (differ) {
    fprintf(stderr, "error:  mer size differ.\n");
    fprintf(stderr, "        AF - "F_U32"\n", ms[0]);
    fprintf(stderr, "        TF - "F_U32"\n", ms[1]);
    fprintf(stderr, "        AC - "F_U32"\n", ms[2]);
    fprintf(stderr, "        DC - "F_U32"\n", ms[3]);
    fprintf(stderr, "        CO - "F_U32"\n", ms[4]);
    exit(1);
  }

  if (dumpFlag) {
    errno = 0;
    dumpSCZF = fopen(dumpSCZFname, "w");
    dumpMCZF = fopen(dumpMCZFname, "w");
    dumpMCSF = fopen(dumpMCSFname, "w");
    dumpMCMF = fopen(dumpMCMFname, "w");
    if (errno)
      fprintf(stderr, "Failed to open the dump files: %s\n", strerror(errno)), exit(1);
  }

  uint32   AFvsAC[NUMCATEGORIES][NUMCATEGORIES];
  uint32   AFvsDC[NUMCATEGORIES][NUMCATEGORIES];
  uint32   AFvsCO[NUMCATEGORIES][NUMCATEGORIES];
  uint32   TFvsAC[NUMCATEGORIES][NUMCATEGORIES];
  uint32   TFvsDC[NUMCATEGORIES][NUMCATEGORIES];
  uint32   TFvsCO[NUMCATEGORIES][NUMCATEGORIES];
  for (uint32 i=0; i<NUMCATEGORIES; i++)
    for (uint32 j=0; j<NUMCATEGORIES; j++) {
      AFvsAC[i][j] = 0;
      AFvsDC[i][j] = 0;
      AFvsCO[i][j] = 0;
      TFvsAC[i][j] = 0;
      TFvsDC[i][j] = 0;
      TFvsCO[i][j] = 0;
    }

  //  The default constructor for kMer sets the mer to size 0, all A.
  //  We need it to be the proper size, and all T.
  kMer   minmer(merSize);

  //  Don't care what we pick, as long as it's a mer in the set.
  //
  if (AF && AF->validMer())  minmer = AF->theFMer();
  if (TF && TF->validMer())  minmer = TF->theFMer();
  if (AC && AC->validMer())  minmer = AC->theFMer();
  if (DC && DC->validMer())  minmer = DC->theFMer();
  if (CO && CO->validMer())  minmer = CO->theFMer();

  speedCounter *C = new speedCounter(" Examining: %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  bool  morestuff = true;
  while (morestuff) {

    //  Find any mer in our set
    if (AF && AF->validMer())  minmer = AF->theFMer();
    if (TF && TF->validMer())  minmer = TF->theFMer();
    if (AC && AC->validMer())  minmer = AC->theFMer();
    if (DC && DC->validMer())  minmer = DC->theFMer();
    if (CO && CO->validMer())  minmer = CO->theFMer();

    //  Find the smallest mer in our set
    if (AF && AF->validMer() && (AF->theFMer() < minmer))  minmer = AF->theFMer();
    if (TF && TF->validMer() && (TF->theFMer() < minmer))  minmer = TF->theFMer();
    if (AC && AC->validMer() && (AC->theFMer() < minmer))  minmer = AC->theFMer();
    if (DC && DC->validMer() && (DC->theFMer() < minmer))  minmer = DC->theFMer();
    if (CO && CO->validMer() && (CO->theFMer() < minmer))  minmer = CO->theFMer();

    //  We need to do up to six comparisons here.
    if (AF && AC)   compare(AF, AC, minmer, AFmode, AFvsAC);
    if (AF && DC)   compare(AF, DC, minmer, AFmode, AFvsDC);
    if (AF && CO)   compare(AF, CO, minmer, AFmode, AFvsCO);
    if (TF && AC)   compare(TF, AC, minmer, TFmode, TFvsAC);
    if (TF && DC)   compare(TF, DC, minmer, TFmode, TFvsDC);
    if (TF && CO)   compare(TF, CO, minmer, TFmode, TFvsCO);

    C->tick();
#if 0
    if (C->tick()) {
      char stringjunk[256];
      fprintf(stderr, "\nMM %s\n", minmer.merToString(stringjunk));
      if (AF) fprintf(stderr, "AF %s\n", AF->theFMer().merToString(stringjunk));
      if (TF) fprintf(stderr, "TF %s\n", TF->theFMer().merToString(stringjunk));
      if (AC) fprintf(stderr, "AC %s\n", AC->theFMer().merToString(stringjunk));
      if (DC) fprintf(stderr, "DC %s\n", DC->theFMer().merToString(stringjunk));
      if (CO) fprintf(stderr, "CO %s\n", CO->theFMer().merToString(stringjunk));
    }
#endif

    //  Advance to the next mer, if we were just used
    morestuff = false;
    if ((AF) && (AF->theFMer() == minmer))   morestuff |= AF->nextMer();
    if ((TF) && (TF->theFMer() == minmer))   morestuff |= TF->nextMer();
    if ((AC) && (AC->theFMer() == minmer))   morestuff |= AC->nextMer();
    if ((DC) && (DC->theFMer() == minmer))   morestuff |= DC->nextMer();
    if ((CO) && (CO->theFMer() == minmer))   morestuff |= CO->nextMer();
  }

  delete C;

  //  output

  if ((AF) && (AC))   output("all frags vs all contigs",          AFmode, AFvsAC);
  if ((AF) && (DC))   output("all frags vs deg. contigs",         AFmode, AFvsDC);
  if ((AF) && (CO))   output("all frags vs non-deg. contigs",     AFmode, AFvsCO);
  if ((TF) && (AC))   output("trimmed frags vs all contigs",      TFmode, TFvsAC);
  if ((TF) && (DC))   output("trimmed frags vs deg. contigs",     TFmode, TFvsDC);
  if ((TF) && (CO))   output("trimmed frags vs non-deg. contigs", TFmode, TFvsCO);

  delete AF;
  delete TF;
  delete AC;
  delete DC;
  delete CO;

  exit(0);
}
