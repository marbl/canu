
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-FEB-20 to 2003-OCT-20
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-30 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-DEC-04 to 2014-APR-11
 *      are Copyright 2005,2007-2008,2011,2013-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "bio++.H"
#include "existDB.H"
#include "libmeryl.H"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

//  Driver for the existDB creation.  Reads a sequence.fasta, builds
//  an existDB for the mers in the file, and then writes the internal
//  structures to disk.
//
//  The existDB constructor is smart enough to read either a pre-built
//  image or a regular multi-fasta file.


int
testFiles(char *filename, char *prefix, uint32 merSize) {
  char *prefixfilename = new char [strlen(prefix) + 32];

  //  Create existDB e and save it to disk
  //
  existDB  *e = new existDB(filename, merSize, existDBnoFlags | existDBcounts, 0, ~uint32ZERO);
  sprintf(prefixfilename, "%s.1", prefix);
  e->saveState(prefixfilename);

  //  Create existDB f by loading the saved copy from disk
  //
  existDB  *f = new existDB(prefixfilename);

  //  Create a fresh existDB g (to check if we corrup the original when saved)
  //
  existDB  *g = new existDB(filename, merSize, existDBnoFlags | existDBcounts, 0, ~uint32ZERO);

  speedCounter *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);
  fprintf(stderr, "Need to iterate over %7.2f Mmers.\n", (uint64MASK(2 * merSize) + 1) / 1000000.0);

  for (uint64 d=0, m=uint64MASK(2 * merSize); m--; ) {
    bool   ee = e->exists(m);
    bool   ef = f->exists(m);
    bool   eg = g->exists(m);

    uint32 ce = e->count(m);
    uint32 cf = f->count(m);
    uint32 cg = g->count(m);

    if ((ee != ef) || (ef != eg) || (ee != eg))
      fprintf(stderr, "mer "uint64HEX" not found : e=%d  f=%d  g=%d\n", m, ee, ef, eg);

    if ((ce != cf) || (cf != cg) || (ce != cg))
      fprintf(stderr, "mer "uint64HEX" count differs : e=%u  f=%u  g=%u (exists=%d)\n", m, ce, cf, cg, ee);

    if ((m & 0xffffff) == 0) {
      //  Been a while since a report, so report.
      d = 1;
    }

    if ((ce > 1) && (d == 1)) {
      //  Report anything not unique, to make sure that we're testing real counts and not just existence.
      fprintf(stderr, "mer "uint64HEX" : e=%u  f=%u  g=%u (exists=%d)\n", m, ce, cf, cg, ee);
      d = 0;
    }

    C->tick();
  }

  delete e;
  delete C;

  return(0);
}


int
testExistence(char *filename, uint32 merSize) {
  existDB         *E      = new existDB(filename, merSize, existDBnoFlags, 0, ~uint32ZERO);
  merStream       *M      = new merStream(new kMerBuilder(merSize), new seqStream(filename), true, true);
  uint64           tried  = 0;
  uint64           lost   = 0;

  while (M->nextMer()) {
    tried++;
    if (!E->exists(M->theFMer()))
      lost++;
  }

  delete M;
  delete E;

  if (lost) {
    fprintf(stderr, "Tried "uint64FMT", didn't find "uint64FMT" merStream mers in the existDB.\n", tried, lost);
    return(1);
  } else {
    return(0);
  }
}



int
testExhaustive(char *filename, char *merylname, uint32 merSize) {
  existDB           *E        = new existDB(filename, merSize, existDBnoFlags, 0, ~uint32ZERO);
  merylStreamReader *M        = new merylStreamReader(merylname);
  speedCounter      *C        = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);
  uint64             found    = uint64ZERO;
  uint64             expected = uint64ZERO;

  FILE              *DUMP     = 0L;

  DUMP = fopen("testExhaustive.ms.dump", "w");

  while (M->nextMer()) {
    if (E->exists(M->theFMer())) {
      expected++;
      fprintf(DUMP, uint64HEX"\n", (uint64)M->theFMer());
    } else {
      fprintf(DUMP, uint64HEX" MISSED!\n", (uint64)M->theFMer());
    }
  }

  fclose(DUMP);

  fprintf(stderr, "Found "uint64FMT" mers in the meryl database.\n", expected);
  fprintf(stderr, "Need to iterate over %7.2f Mmers.\n", (uint64MASK(2 * merSize) + 1) / 1000000.0);

  DUMP = fopen("testExhaustive.ck.dump", "w");

  for (uint64 m = uint64MASK(2 * merSize); m--; ) {
    if (E->exists(m)) {
      found++;
      fprintf(DUMP, uint64HEX"\n", m);
    }
    C->tick();
  }

  fclose(DUMP);

  delete C;
  delete E;
  delete M;

  if (expected != found) {
    fprintf(stderr, "Expected to find "uint64FMT" mers, but found "uint64FMT" instead.\n",
            expected, found);
    return(1);
  } else {
    return(0);
  }
}


const char *usage =
"usage: %s [stuff]\n"
"       -mersize mersize\n"
"         -- Use the specified mersize when building existDB tables.\n"
"\n"
"       -build some.fasta prefix\n"
"         -- Build an existDB on all mers in some.fasta and save\n"
"            the tables into prefix.\n"
"\n"
"       -describe prefix\n"
"         -- Reports the state of some existDB file.\n"
"\n"
"       -testfiles some.fasta prefix\n"
"         -- Build an existDB table from some.fasta.  Write that table to disk.\n"
"            Load the table back.  Compare that each mer in some.fasta is present\n"
"            in all three existDB tables created earlier.\n"
"\n"
"       -testexistence some.fasta\n"
"         -- Build an existDB table from some.fasta, check that every\n"
"            mer in some.fasta can be found in the table.  Does not\n"
"            guarantee that every mer in the table is found in the file.\n"
"\n"
"       -testexhaustive some.fasta some.meryl\n"
"         -- Build an existDB table from some.fasta, check _EVERY_ mer\n"
"            for existance.  Complain if a mer exists in the table but\n"
"            not in the meryl database.  Assumes 'some.meryl' is the\n"
"            mercount of some.fasta.\n"
"\n";

int
main(int argc, char **argv) {
  uint32    mersize = 20;

  if (argc < 3) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-mersize", 2) == 0) {
      arg++;
      mersize = atoi(argv[arg]);

    } else if (strncmp(argv[arg], "-describe", 2) == 0) {
      existDB *e = new existDB(argv[argc-1], false);
      e->printState(stdout);
      delete e;
      exit(0);

    } else if (strncmp(argv[arg], "-testfiles", 8) == 0) {
      exit(testFiles(argv[arg+1], argv[arg+2], mersize));

    } else if (strncmp(argv[arg], "-testexistence", 8) == 0) {
      exit(testExistence(argv[arg+1], mersize));

    } else if (strncmp(argv[arg], "-testexhaustive", 8) == 0) {
      exit(testExhaustive(argv[arg+1], argv[arg+2], mersize));

    } else if (strncmp(argv[arg], "-build", 2) == 0) {
      existDB  *e = new existDB(argv[argc-2], mersize, existDBnoFlags, 0, ~uint32ZERO);
      e->saveState(argv[argc-1]);
      delete e;
      exit(0);
    }

    arg++;
  }

  exit(0);
}
