#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bio++.H"

//  Build a merStreamFile using small mers, read it back using bigger mers.

//  Also report the correct stuff.
//#define SHOW_CORRECT

#define BUILD_SIZE   14
#define TEST_SIZE    31

#if KMER_WORDS * 32 < TEST_SIZE
#error KMER_WORDS too small for TEST_SIZE
#endif



int
test(u64bit compared, merStream *M, merStreamFileReader *R) {
  static char  stra[256], strb[256], strc[256], strd[256];

  if ((M->theFMer()                != R->theFMer()) ||
      (M->theRMer()                != R->theRMer()) ||
      (M->thePositionInSequence()  != R->thePositionInSequence()) ||
      (M->thePositionInStream()    != R->thePositionInStream()) ||
      (M->theSequenceNumber()      != R->theSequenceNumber())) {
    fprintf(stderr, u64bitFMTW(8)": !!! M got F=%s R=%s "u64bitFMTW(7)" "u64bitFMTW(7)" "u64bitFMTW(7)" but R got F=%s R=%s "u64bitFMTW(7)" "u64bitFMTW(7)" "u64bitFMTW(7)"\n",
            compared,
            M->theFMer().merToString(stra), M->theRMer().merToString(strb), M->thePositionInSequence(), M->thePositionInStream(), M->theSequenceNumber(),
            R->theFMer().merToString(strc), R->theRMer().merToString(strd), R->thePositionInSequence(), R->thePositionInStream(), R->theSequenceNumber());
    
    return(1);
  }

#ifdef SHOW_CORRECT
  fprintf(stderr, u64bitFMTW(8)":     M got F=%s R=%s "u64bitFMTW(7)" "u64bitFMTW(7)" "u64bitFMTW(7)" and R got F=%s R=%s "u64bitFMTW(7)" "u64bitFMTW(7)" "u64bitFMTW(7)"\n",
          compared,
          M->theFMer().merToString(stra), M->theRMer().merToString(strb), M->thePositionInSequence(), M->thePositionInStream(), M->theSequenceNumber(),
          R->theFMer().merToString(strc), R->theRMer().merToString(strd), R->thePositionInSequence(), R->thePositionInStream(), R->theSequenceNumber());
#endif

  return(0);
}



int
streamingTest(char *msfile, bool inCore) {

  fprintf(stderr, "streamingTest-- %s %s\n", msfile, (inCore) ? "inCore" : "disk-based");

  int                  errors = 0;
  u64bit               compared = 0;
  merStreamFileReader *R = new merStreamFileReader("junk", TEST_SIZE);
  FastAstream         *F = new FastAstream(msfile);
  merStream           *M = new merStream(TEST_SIZE, F);

  if (inCore)
    R->loadInCore();

  double startTime = getTime();

  fprintf(stderr, "Testing streaming access.\n");
  while (M->nextMer() && R->nextMer()) {
    compared++;
    errors += test(compared, M, R);
  }

  fprintf(stderr, "Compared "u64bitFMT" mers. (%f seconds)\n", compared, getTime() - startTime);

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
randomAccessTest(char *msfile, bool inCore) {

  fprintf(stderr, "randomAccessTest -- %s %s\n", msfile, (inCore) ? "inCore" : "disk-based");

  int                  errors = 0;
  u64bit               merNum   = 0;
  merStreamFileReader *R = new merStreamFileReader("junk", TEST_SIZE);
  FastAstream         *F = new FastAstream(msfile);
  merStream           *M = new merStream(TEST_SIZE, F);

  if (inCore)
    R->loadInCore();

  u64bit numMers = R->numberOfMers();

  double startTime = 0.0;
  double seekTime  = 0.0;

  //  Load the first mer from the stream.
  M->nextMer();

  //  How many seeks?  At most 100, but we'll decrease if the seeks
  //  are less than 1000 bases.
  //
  u64bit numSeeks = 100000;
  while (numMers / numSeeks < 1000)
    numSeeks = (u64bit)floor(numSeeks * 0.8);

  for (u64bit s=numSeeks; --s; ) {

    //  Skip 'skipSize' mers in M
    //
    for (u32bit loop=0; loop < numMers / numSeeks; loop++, merNum++)
      M->nextMer();

    //  Seek to the proper spot in R
    //
    startTime = getTime();
    R->setIterationStart(merNum);
    R->nextMer();
    seekTime += getTime() - startTime;

    errors += test(merNum, M, R);
  }

  fprintf(stderr, "(%f seconds)\n", seekTime);

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

  fprintf(stderr, "Building merStreamFile.\n");

  merStreamFileBuilder   *B = new merStreamFileBuilder(BUILD_SIZE, argv[1], "junk");
  u64bit numMers = B->build(true);
  delete B;

  fprintf(stderr, "Found "u64bitFMT" mers in %s at mersize %d\n", numMers, argv[1], BUILD_SIZE);

  u32bit errors = 0;
  errors += streamingTest(argv[1], false);
  errors += streamingTest(argv[1], true);
  errors += randomAccessTest(argv[1], false);
  errors += randomAccessTest(argv[1], true);

  if (errors > 0) {
    fprintf(stderr, "\nThere were "u32bitFMT" errors.\n", errors);
    exit(1);
  }

  unlink("junk.merStream");

  exit(0);
}
