#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include "positionDB.H"
#include "fasta.H"


#define MERSIZE   20
#define TABLESIZE 20


char *usage = 
"usage: %s -d <a.fasta>\n"
"usage: %s -a <a.fasta>\n"
"usage: %s -s <1.fasta> <2.fasta>\n"
"\n"
"-d:  Builds a positionDB and dumps it.\n"
"\n"
"-a:  Tests if each and every mer is found in the positionDB.\n"
"     Reports if it doesn't find a mer at the correct position.\n"
"     Doesn't report if table has too much stuff.\n"
"\n"
"-s:  Builds a positionDB from 1.fasta, does lookup of 2.fasta.\n"
"     Reports if it doesn't find a mer.\n";


void
buildOnly(char *filename) {
  FastA       F(filename, false, false);
  FastABuffer B;

  F.first(B);

  positionDB *M = new positionDB(B.sequence(), MERSIZE, TABLESIZE, true);
  //M->dumpTable();
  delete M;
}


void
dump(char *filename) {
  FastA       F(filename, false, false);
  FastABuffer B;

  F.first(B);

  positionDB *M = new positionDB(B.sequence(), MERSIZE, TABLESIZE, true);
  M->dumpTable();
  delete M;
}


void
all(char *filename) {
  FastA        O(filename, false, false);
  FastABuffer  B;

  O.first(B);

  positionDB *DB = new positionDB(B.sequence(), MERSIZE, TABLESIZE, true);
  u64bit     *posn    = new u64bit [1024];
  u64bit      posnMax = 1024;
  u64bit      posnLen = 0;

  merStream *T = new merStream(MERSIZE, filename);

  for (u32bit j=0; T->nextMer(); j++) {
    if (DB->get(T->theFMer(),
                posn,
                posnMax,
                posnLen)) {
      bool missing = true;
      for (u32bit i=0; i<posnLen; i++)
        if (posn[i] == T->thePosition())
          missing = false;
      if (missing) {
        fprintf(stdout, "Found table entries (%lu), but no position for mer=%s at pos=%d\n",
                posnLen, T->theFMerString(), T->thePosition());
        fflush(stdout);
#ifdef ERROR
        fprintf(stderr, "Found table entries (%lu), but no position for mer=%s at pos=%d\n",
                posnLen, T->theFMerString(), T->thePosition());
        fflush(stderr);
#endif
      }
    } else {
      fprintf(stdout, "Found no matches for mer=%s at pos=%d\n",
              T->theFMerString(), T->thePosition());
      fflush(stdout);
#ifdef ERROR
      fprintf(stderr, "Found no matches for mer=%s at pos=%d\n",
              T->theFMerString(), T->thePosition());
      fflush(stderr);
#endif
    }
  }
}

void
search(char *filename, char *query) {
  FastA        O(filename, false, false);
  FastABuffer  B;

  O.first(B);

  positionDB *DB = new positionDB(B.sequence(), MERSIZE, TABLESIZE, true);
  u64bit     *posn    = new u64bit [1024];
  u64bit      posnMax = 1024;
  u64bit      posnLen = 0;

  merStream *T = new merStream(MERSIZE, query);

  for (u32bit j=0; T->nextMer(); j++) {
    if (DB->get(T->theFMer(),
                posn,
                posnMax,
                posnLen)) {
      fprintf(stdout, "Got a F match for mer=%s at onePos=%d (in mers), numMatches=%ld\n",
              T->theFMerString(), j, posnLen);
      fflush(stdout);
#ifdef ERROR
      fprintf(stderr, "Got a F match for mer=%s at onePos=%d (in mers), numMatches=%ld\n",
              T->theFMerString(), j, posnLen);
      fflush(stderr);
#endif
    }

    if (DB->get(T->theRMer(),
                posn,
                posnMax,
                posnLen)) {
      fprintf(stdout, "Got a R match for mer=%s at onePos=%d (in mers), numMatches=%ld\n",
              T->theRMerString(), j, posnLen);
      fflush(stdout);
#ifdef ERROR
      fprintf(stderr, "Got a R match for mer=%s at onePos=%d (in mers), numMatches=%ld\n",
              T->theRMerString(), j, posnLen);
      fflush(stderr);
#endif
    }
  }
}






int
main(int argc, char **argv) {

  if (argc == 1) {
    fprintf(stderr, usage, argv[0], argv[0], argv[0]);
    exit(1);
  }

  if (argv[1][0] == '-') {
    switch (argv[1][1]) {
      case 'b':
        buildOnly(argv[2]);
        break;
      case 'd':
        dump(argv[2]);
        break;
      case 'a':
        all(argv[2]);
        break;
      case 's':
        search(argv[2], argv[3]);
        break;
    }
  } else {
    fprintf(stderr, usage, argv[0], argv[0], argv[0]);
    exit(1);
  }
}

