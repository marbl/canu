#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "bio++.H"
#include "positionDB.H"

//  Driver for the positionDB creation.  Reads a sequence.fasta, builds
//  a positionDB for the mers in the file, and then writes the internal
//  structures to disk.
//
//  The positionDB constructor is smart enough to read either a pre-built
//  image or a regular multi-fasta file.


#define MERSIZE 20
#define TBLSIZE 19

char const *usage =
"usage: %s [-mersize s] [-merskip k] [-tablesize t]\n"
"         default: s=20, k=0, t=auto\n"
"\n"
"       The default operation is to create an image:\n"
"         [build opts] <sequence.fasta> <datafile>\n"
"\n"
"       To dump information about an image:\n"
"         -dump datafile\n"
"\n"
"       To run sanity tests:\n"
"         -buildonly [build opts] sequence.fasta\n"
"           --  just builds a table and exits\n"
"\n"
"         -existence [build opts] sequence.fasta\n"
"           --  builds (or reads) a table reports if any mers\n"
"               in sequence.fasta cannot be found\n"
"\n"
"         -extra [build opts] sequence.fasta\n"
"           --  builds (or reads) a table reports if any mers\n"
"               NOT in sequence.fasta are be found\n"
"\n"
"         -test1 sequence.fasta\n"
"           --  Tests if each and every mer is found in the\n"
"               positionDB.  Reports if it doesn't find a mer\n"
"               at the correct position.  Doesn't report if table\n"
"               has too much stuff.\n"
"\n"
"         -test2 db.fasta sequence.fasta\n"
"           --  Builds a positionDB from db.fasta, then searches\n"
"               the table for each mer in sequence.fasta.  Reports\n"
"               all mers it finds.\n"
"            -- This is a silly test and you shouldn't do it.\n";



int
test1(char *filename) {
  FastAstream *F       = new FastAstream(filename);
  merStream   *T       = new merStream(MERSIZE, F);
  positionDB  *M       = new positionDB(T, MERSIZE, 0, TBLSIZE, 0L, 0L, true);
  u64bit      *posn    = new u64bit [1024];
  u64bit       posnMax = 1024;
  u64bit       posnLen = u64bitZERO;
  u32bit       missing = u32bitZERO;
  u32bit       failed  = u32bitZERO;
  char         str[33];

  if (T->rewind() == false) {
    fprintf(stderr, "test 1 failed to rewind the merStream.\n");
    return(true);
  }

  while (T->nextMer()) {
    if (M->get(T->theFMer(),
                posn,
                posnMax,
                posnLen)) {

      missing = u32bitZERO;
      for (u32bit i=0; i<posnLen; i++)
        if (posn[i] == T->thePositionInStream())
          missing++;

      if (missing != 1) {
        failed++;

        fprintf(stdout, "%s @ "u64bitFMT"/"u64bitFMT": Found "u64bitFMT" table entries, and "u32bitFMT" matching positions (",
                T->theFMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen, missing);

        for (u32bit i=0; i<posnLen; i++) {
          fprintf(stdout, u64bitFMT, posn[i]);
          if (i < posnLen - 1)
            fprintf(stdout, " ");
          else
            fprintf(stdout, ")\n");
        }
      }
    } else {
      failed++;

      fprintf(stdout, "Found no matches for mer=%s at pos="u64bitFMT"\n",
              T->theFMer().merToString(str), T->thePositionInStream());
    }
  }

  delete M;
  delete T;
  delete F;

  return(failed != 0);
}



int
test2(char *filename, char *query) {
  FastAstream *F       = new FastAstream(filename);
  merStream   *T       = new merStream(MERSIZE, F);
  positionDB  *M       = new positionDB(T, MERSIZE, 0, TBLSIZE, 0L, 0L, true);
  u64bit      *posn    = new u64bit [1024];
  u64bit       posnMax = 1024;
  u64bit       posnLen = u64bitZERO;
  char         str[33];

  delete T;
  delete F;

  F = new FastAstream(query);
  T = new merStream(MERSIZE, F);

  while (T->nextMer()) {
    if (M->get(T->theFMer(),
                posn,
                posnMax,
                posnLen)) {
      fprintf(stdout, "Got a F match for mer=%s at "u64bitFMT"/"u64bitFMT" (in mers), numMatches="u64bitFMT"\n",
              T->theFMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen);
    }

    if (M->get(T->theRMer(),
                posn,
                posnMax,
                posnLen)) {
      fprintf(stdout, "Got a R match for mer=%s at "u64bitFMT"/"u64bitFMT" (in mers), numMatches="u64bitFMT"\n",
              T->theRMer().merToString(str), T->theSequenceNumber(), T->thePositionInStream(), posnLen);
    }
  }

  delete F;
  delete M;
  delete T;

  return(0);
}



int
main(int argc, char **argv) {
  u32bit    mersize = 20;
  u32bit    merskip = 0;
  u32bit    tblsize = 0;

  if (argc < 3) {
    fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-mersize", 6) == 0) {
      arg++;
      mersize = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-merskip", 6) == 0) {
      arg++;
      merskip = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-tablesize", 3) == 0) {
      arg++;
      tblsize = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-dump", 2) == 0) {
      positionDB *e = new positionDB(argv[argc-1], false);
      e->printState(stdout);
      delete e;
      exit(0);
    } else if (strncmp(argv[arg], "-test1", 6) == 0) {
      exit(test1(argv[arg+1]));
    } else if (strncmp(argv[arg], "-test2", 6) == 0) {
      exit(test2(argv[arg+1], argv[arg+2]));
    }
    arg++;
  }

  fprintf(stderr, "Building table with merSize "u32bitFMT", merSkip "u32bitFMT" and table size "u32bitFMT"\n", mersize, merskip, tblsize);

  u64bit    seqlen  = 0;
  u32bit    ts      = 25;

  fprintf(stderr, "Opening '%s', finding length of sequences\n", argv[argc-2]);
  FastAWrapper  *F = new FastAWrapper(argv[argc-2]);
  F->openIndex();
  for (u32bit i=0; i<F->getNumberOfSequences(); i++)
    seqlen += F->sequenceLength(i);

  if (seqlen < 64 * 1024 * 1024) ts = 24;
  if (seqlen < 16 * 1024 * 1024) ts = 23;
  if (seqlen <  4 * 1024 * 1024) ts = 22;
  if (seqlen <  2 * 1024 * 1024) ts = 21;
  if (seqlen <  1 * 1024 * 1024) ts = 20;

  if (tblsize == 0)
    tblsize = ts;

  if (ts != tblsize)
    fprintf(stderr, "For sequence length "u64bitFMT", suggested tableSize is "u32bitFMT", you want "u32bitFMT"\n", seqlen, ts, tblsize);

#if 0
  positionDB  *e = new positionDB(0L, argv[argc-2], mersize, merskip, tblsize, true);

  fprintf(stderr, "Writing positionDB state to '%s'\n", argv[argc-1]);
  e->saveState(argv[argc-1]);

  delete e;
#endif

  exit(0);
}
