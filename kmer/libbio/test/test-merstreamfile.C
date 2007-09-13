#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bio++.H"

//  Build a merStreamFile using small mers, read it back using bigger mers.

//  Also report the correct stuff.
#undef SHOW_CORRECT

#define BUILD_SIZE   14
#define TEST_SIZE    31

#if KMER_WORDS * 32 < TEST_SIZE
#error KMER_WORDS too small for TEST_SIZE
#endif



int
test(u64bit compared, merStream *M, merStream *R) {
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

  merStream    *merstr = new merStream(new kMerBuilder(TEST_SIZE), new seqStream(msfile, true));
  seqStore     *seqsto = new seqStore(msfile, new seqStream(msfile, true));

  if (inCore)
    seqsto->loadStoreInCore();

  if (0) {
    char stra[256];

    while (merstr->nextMer()) {
      fprintf(stderr, "seqpos:"u64bitFMT" seqiid:"u64bitFMT" strpos:"u64bitFMT" %s\n",
              merstr->thePositionInSequence(), merstr->theSequenceNumber(), merstr->thePositionInStream(),
              merstr->theFMer().merToString(stra));
    }
    delete merstr;
    merstr = new merStream(new kMerBuilder(TEST_SIZE), new seqStream(msfile, true));
    exit(1);
  }

  //  Just dump the sequence
  if (0) {
    char ch;
    char st[1024] = {0};
    int  i = 0;
    while (ch = seqsto->get()) {
      st[i++] = ch;
    }
    printf("%s\n", st);

    seqsto->rewind();
  }

  merStream    *mersto = new merStream(new kMerBuilder(TEST_SIZE), seqsto);

  double startTime = getTime();

  fprintf(stderr, "Testing streaming access.\n");
  while (merstr->nextMer() && mersto->nextMer()) {
    compared++;
    errors += test(compared, merstr, mersto);
  }

  fprintf(stderr, "Compared "u64bitFMT" mers. (%f seconds)\n", compared, getTime() - startTime);

  if (merstr->nextMer()) {
    fprintf(stderr, "ERROR: Extra mers in the merStream!\n");
    errors++;
  }

  if (mersto->nextMer()) {
    fprintf(stderr, "ERROR: Extra mers in the merStore!\n");
    errors++;
  }

  delete mersto;
  delete merstr;

  return(errors);
}



int
randomAccessTest(char *msfile, bool inCore) {
  int                  errors = 0;

  fprintf(stderr, "randomAccessTest skipped until it gets rewritten.\n");

#if 0
  fprintf(stderr, "randomAccessTest -- %s %s\n", msfile, (inCore) ? "inCore" : "disk-based");

  u64bit               merNum   = 0;

  merStream           *merstr = new merStream(new kMerBuilder(TEST_SIZE), new seqStream(msfile, true));
  seqStore            *seqsto = new seqStore(msfile, new seqStream(msfile, true));
  merStream           *mersto = new merStream(new kMerBuilder(TEST_SIZE), seqsto);

  if (inCore)
    seqsto->loadStoreInCore();

  u64bit numMers = seqsto->numberOfACGT();

  double startTime = 0.0;
  double seekTime  = 0.0;

  //  Load the first mer from the stream.
  merstr->nextMer();

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
      merstr->nextMer();

    //  Seek to the proper spot in R
    //
    startTime = getTime();
    mersto->seek(merNum);
    seekTime += getTime() - startTime;

    errors += test(merNum, merstr, mersto);

    if (errors > 20)
      exit(1);
  }

  fprintf(stderr, "(%f seconds)\n", seekTime);

  delete mersto;
  delete merstr;
  delete seqsto;
#endif

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

  seqStream *seqstr = new seqStream(argv[1], true);
  seqStore  *seqsto = new seqStore(argv[1], seqstr);

  fprintf(stderr, "Found "u64bitFMT" ACGT in %s\n", seqsto->numberOfACGT(), argv[1]);

  delete seqsto;
  delete seqstr;

  u32bit errors = 0;
  errors += streamingTest(argv[1], false);
  errors += streamingTest(argv[1], true);
  errors += randomAccessTest(argv[1], false);
  errors += randomAccessTest(argv[1], true);

  if (errors > 0) {
    fprintf(stderr, "\nThere were "u32bitFMT" errors.\n", errors);
    exit(1);
  }

  unlink("junk.seqStream.blocks");
  unlink("junk.seqStream.sequence");

  exit(0);
}
