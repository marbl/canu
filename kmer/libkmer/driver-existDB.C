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
testFiles(char *filename, char *prefix, u32bit merSize) {
  char *prefixfilename = new char [strlen(prefix) + 32];

  //  Create existDB e and save it to disk
  //
  existDB  *e = new existDB(filename, merSize, existDBnoFlags, 0, ~u32bitZERO);
  sprintf(prefixfilename, "%s.1", prefix);
  e->saveState(prefixfilename);

  //  Create existDB f by loading the saved copy from disk
  //
  existDB  *f = new existDB(prefixfilename);

  //  Create a fresh existDB g
  //
  existDB  *g = new existDB(filename, merSize, existDBnoFlags, 0, ~u32bitZERO);

  speedCounter *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);
  fprintf(stderr, "Need to iterate over %7.2f Mmers.\n", (u64bitMASK(2 * merSize) + 1) / 1000000.0);

  for (u64bit m = u64bitMASK(2 * merSize); m--; ) {
    bool ee = e->exists(m);
    bool ef = f->exists(m);
    bool eg = g->exists(m);

    if ((ee != ef) || (ef != eg) || (ee != eg))
      fprintf(stderr, "mer "u64bitHEX" not found : e=%d  f=%d  g=%d\n", m, ee, ef, eg);

    C->tick();
  }

  delete e;
  delete C;

  return(0);
}


int
testExistence(char *filename, u32bit merSize) {
  existDB         *E      = new existDB(filename, merSize, existDBnoFlags, 0, ~u32bitZERO);
  merStream       *M      = new merStream(new kMerBuilder(merSize), new seqStream(filename), true, true);
  u64bit           tried  = 0;
  u64bit           lost   = 0;

  while (M->nextMer()) {
    tried++;
    if (!E->exists(M->theFMer()))
      lost++;
  }

  delete M;
  delete E;

  if (lost) {
    fprintf(stderr, "Tried "u64bitFMT", didn't find "u64bitFMT" merStream mers in the existDB.\n", tried, lost);
    return(1);
  } else {
    return(0);
  }
}



int
testExhaustive(char *filename, char *merylname, u32bit merSize) {
  existDB           *E        = new existDB(filename, merSize, existDBnoFlags, 0, ~u32bitZERO);
  merylStreamReader *M        = new merylStreamReader(merylname);
  speedCounter      *C        = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);
  u64bit             found    = u64bitZERO;
  u64bit             expected = u64bitZERO;

  FILE              *DUMP     = 0L;

  DUMP = fopen("testExhaustive.ms.dump", "w");

  while (M->nextMer()) {
    if (E->exists(M->theFMer())) {
      expected++;
      fprintf(DUMP, "0x%016lx\n", (u64bit)M->theFMer());
    } else {
      fprintf(DUMP, "0x%016lx MISSED!\n", (u64bit)M->theFMer());
    }
  }

  fclose(DUMP);

  fprintf(stderr, "Found "u64bitFMT" mers in the meryl database.\n", expected);
  fprintf(stderr, "Need to iterate over %7.2f Mmers.\n", (u64bitMASK(2 * merSize) + 1) / 1000000.0);

  DUMP = fopen("testExhaustive.ck.dump", "w");

  for (u64bit m = u64bitMASK(2 * merSize); m--; ) {
    if (E->exists(m)) {
      found++;
      fprintf(DUMP, "0x%016lx\n", m);
    }
    C->tick();
  }

  fclose(DUMP);

  delete C;
  delete E;
  delete M;

  if (expected != found) {
    fprintf(stderr, "Expected to find "u64bitFMT" mers, but found "u64bitFMT" instead.\n",
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
  u32bit    mersize = 20;

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
      existDB  *e = new existDB(argv[argc-2], mersize, existDBnoFlags, 0, ~u32bitZERO);
      e->saveState(argv[argc-1]);
      delete e;
      exit(0);
    }
    arg++;
  }

  exit(0);
}
