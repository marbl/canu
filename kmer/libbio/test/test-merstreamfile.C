#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bio++.H"

int
streamingTest(char *msfile) {
  int                  errors = 0;
  u64bit               compared = 0;
  merStreamFileReader *R = new merStreamFileReader("junk");
  FastAstream         *F = new FastAstream(msfile);
  merStream           *M = new merStream(20, F);

  fprintf(stderr, "Testing streaming access.\n");
  while (M->nextMer() && R->nextMer()) {
    compared++;
    if ((M->theFMer()                != R->theFMer()) ||
        (M->theRMer()                != R->theRMer()) ||
        (M->thePositionInSequence()  != R->thePositionInSequence()) ||
        (M->thePositionInStream()    != R->thePositionInStream()) ||
        (M->theSequenceNumber()      != R->theSequenceNumber())) {
      fprintf(stderr, u64bitFMT": !!! M got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT" but R got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT"\n",
              compared,
              M->theFMer(), M->theRMer(), M->thePositionInSequence(), M->thePositionInStream(), M->theSequenceNumber(),
              R->theFMer(), R->theRMer(), R->thePositionInSequence(), R->thePositionInStream(), R->theSequenceNumber());
      errors++;
#if DEBUG
    } else {
      fprintf(stderr, u64bitFMT":     M got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT" but R got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT"\n",
              compared,
              M->theFMer(), M->theRMer(), M->thePositionInSequence(), M->thePositionInStream(), M->theSequenceNumber(),
              R->theFMer(), R->theRMer(), R->thePositionInSequence(), R->thePositionInStream(), R->theSequenceNumber());
#endif
    }
  }

  fprintf(stderr, "Compared "u64bitFMT" mers.\n", compared);

  if (M->nextMer()) {
    fprintf(stderr, "ERROR: Extra mers in the merStream!\n");
    errors++;
  }

  if (R->nextMer()) {
    fprintf(stderr, "ERROR: Extra mers in the merStreamFile!\n");
    errors++;
  }

  delete M;
  delete F;
  delete R;

  return(errors);
}



int
randomAccessTest(char *msfile, u64bit numMers) {
  int                  errors = 0;
  u64bit               merNum   = 0;
  merStreamFileReader *R = new merStreamFileReader("junk");
  FastAstream         *F = new FastAstream(msfile);
  merStream           *M = new merStream(20, F);

  //  Load the first mer from the stream.
  M->nextMer();

  //  How many seeks?  At most 100, but we'll decrease if the seeks
  //  are less than 1000 bases.
  //
  u64bit numSeeks = 100;
  while (numMers / numSeeks < 1000)
    numSeeks = (u64bit)floor(numSeeks * 0.8);

  fprintf(stderr, "Testing random access on "u64bitFMT" seeks of size "u64bitFMT".\n", numSeeks, numMers / numSeeks);
  for (u64bit s=numSeeks; --s; ) {
    fprintf(stderr, " "u64bitFMT" seeks remain\r", s);
    fflush(stderr);

    //  Skip 'skipSize' mers in M
    //
    for (u32bit loop=0; loop < numMers / numSeeks; loop++, merNum++)
      M->nextMer();

    //  Seek to the proper spot in R
    //
    R->seekToMer(merNum);
    R->nextMer();

    if ((M->theFMer()               != R->theFMer()) ||
        (M->theRMer()               != R->theRMer()) ||
        (M->thePositionInSequence() != R->thePositionInSequence()) ||
        (M->thePositionInStream()   != R->thePositionInStream()) ||
        (M->theSequenceNumber()     != R->theSequenceNumber())) {
      fprintf(stderr, u64bitFMT": !!! M got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT" but R got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT"\n",
              merNum,
              M->theFMer(), M->theRMer(), M->thePositionInSequence(), M->thePositionInStream(), M->theSequenceNumber(),
              R->theFMer(), R->theRMer(), R->thePositionInSequence(), R->thePositionInStream(), R->theSequenceNumber());
      errors++;
#if DEBUG
    } else {
      fprintf(stderr, u64bitFMT":     M got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT" but R got F="u64bitHEX" R="u64bitHEX" "u64bitFMT" "u64bitFMT" "u64bitFMT"\n",
              merNum,
              M->theFMer(), M->theRMer(), M->thePositionInSequence(), M->thePositionInStream(), M->theSequenceNumber(),
              R->theFMer(), R->theRMer(), R->thePositionInSequence(), R->thePositionInStream(), R->theSequenceNumber());
#endif
    }
  }

  delete M;
  delete F;
  delete R;

  return(errors);
}



int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "       Builds a merStreamFile, and then checks that it returns\n");
    fprintf(stderr, "       exactly the same stuff as a merStream(\"some.fasta\") does.\n");
    fprintf(stderr, "       Returns 1 if error, 0 if OK\n");
    exit(1);
  }

  //  Build a compressed merstreamfile.
  //
  merStreamFileBuilder   *B = new merStreamFileBuilder(20, argv[1], "junk");
  u64bit numMers = B->build(true);
  delete B;

  fprintf(stderr, "Found "u64bitFMT" mers in %s\n", numMers, argv[1]);

  u32bit errors = 0;
  errors += streamingTest(argv[1]);
  errors += randomAccessTest(argv[1], numMers);

  if (errors > 0) {
    fprintf(stderr, "\nThere were "u32bitFMT" errors.\n", errors);
    exit(1);
  }

  unlink("junk.merStream");

  exit(0);
}
