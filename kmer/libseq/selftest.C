#include "util.h"

#include "seqFile.H"
#include "sffFile.H"
#include "fastaFile.H"
#include "seqStore.H"
//#include "gkpStore.H"
#include "seqFactory.H"
#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

//  Include all header files directly; this is test code, right?


void
testSSsimple(seqStream *SS) {

  fprintf(stdout, "testSSsimple() begins.\n");

  while (!SS->eof()) {
    u64bit  seqpos = SS->seqPos();
    u64bit  seqiid = SS->seqIID();
    u64bit  strpos = SS->strPos();
    char    x      = SS->get();

    fprintf(stdout, "%c seqPos="u64bitFMT" seqIID="u64bitFMT" strPos="u64bitFMT"\n",
            (x) ? x : '0', seqpos, seqiid, strpos);
  }

  fprintf(stdout, "testSSsimple() ends.\n");
}


void
testMSsimple(merStream *MS) {
  char  testmer[32];
  bool  verbose = false;
  bool  nm      = false;

  if (verbose)
    fprintf(stdout, "testMSsimple() begins.\n");

  nm = MS->nextMer();
  assert(nm == true);
  assert(MS->thePositionInSequence() == 0);
  assert(MS->thePositionInStream()   == 0);
  assert(MS->theSequenceNumber()     == 0);
  assert(strcmp(MS->theFMer().merToString(testmer), "GGGTCAACTCCGCCCGCACT") == 0);

  if (verbose) {
    fprintf(stdout, "MS 1 posInSeq="u64bitFMT" posInStr="u64bitFMT" seqNum="u64bitFMT"\n",
            MS->thePositionInSequence(),
            MS->thePositionInStream(),
            MS->theSequenceNumber());
    if (strcmp(MS->theFMer().merToString(testmer), "GGGTCAACTCCGCCCGCACT"))
      fprintf(stdout, "MS 1 failed '%s' != ''.\n", testmer);
  }

  nm = MS->nextMer();
  assert(nm == true);
  assert(MS->thePositionInSequence() == 1);
  assert(MS->thePositionInStream()   == 1);
  assert(MS->theSequenceNumber()     == 0);
  assert(strcmp(MS->theFMer().merToString(testmer), "GGTCAACTCCGCCCGCACTC") == 0);

  if (verbose) {
    fprintf(stdout, "MS 2 posInSeq="u64bitFMT" posInStr="u64bitFMT" seqNum="u64bitFMT"\n",
            MS->thePositionInSequence(),
            MS->thePositionInStream(),
            MS->theSequenceNumber());
    if (strcmp(MS->theFMer().merToString(testmer), "GGTCAACTCCGCCCGCACTC"))
      fprintf(stdout, "MS 2 failed '%s' != ''.\n", testmer);
  }

  nm = MS->nextMer();
  assert(nm == true);
  assert(MS->thePositionInSequence() == 2);
  assert(MS->thePositionInStream()   == 2);
  assert(MS->theSequenceNumber()     == 0);
  assert(strcmp(MS->theFMer().merToString(testmer), "GTCAACTCCGCCCGCACTCT") == 0);

  if (verbose) {
    fprintf(stdout, "MS 3 posInSeq="u64bitFMT" posInStr="u64bitFMT" seqNum="u64bitFMT"\n",
            MS->thePositionInSequence(),
            MS->thePositionInStream(),
            MS->theSequenceNumber());
    if (strcmp(MS->theFMer().merToString(testmer), "GTCAACTCCGCCCGCACTCT"))
      fprintf(stdout, "MS 3 failed '%s' != ''.\n", testmer);
  }

  nm = MS->nextMer();
  assert(nm == true);
  assert(MS->thePositionInSequence() == 3);
  assert(MS->thePositionInStream()   == 3);
  assert(MS->theSequenceNumber()     == 0);
  assert(strcmp(MS->theFMer().merToString(testmer), "TCAACTCCGCCCGCACTCTA") == 0);

  if (verbose) {
    fprintf(stdout, "MS 4 posInSeq="u64bitFMT" posInStr="u64bitFMT" seqNum="u64bitFMT"\n",
            MS->thePositionInSequence(),
            MS->thePositionInStream(),
            MS->theSequenceNumber());
    if (strcmp(MS->theFMer().merToString(testmer), "TCAACTCCGCCCGCACTCTA"))
      fprintf(stdout, "MS 4 failed '%s' != ''.\n", testmer);
  }

  nm = MS->nextMer();
  assert(nm == true);
  assert(MS->thePositionInSequence() == 4);
  assert(MS->thePositionInStream()   == 4);
  assert(MS->theSequenceNumber()     == 0);
  assert(strcmp(MS->theFMer().merToString(testmer), "CAACTCCGCCCGCACTCTAG") == 0);

  if (verbose) {
    fprintf(stdout, "MS 5 posInSeq="u64bitFMT" posInStr="u64bitFMT" seqNum="u64bitFMT"\n",
            MS->thePositionInSequence(),
            MS->thePositionInStream(),
            MS->theSequenceNumber());
    if (strcmp(MS->theFMer().merToString(testmer), "CAACTCCGCCCGCACTCTAG"))
      fprintf(stdout, "MS 5 failed '%s' != ''.\n", testmer);
  }

  nm = MS->nextMer();
  assert(nm == true);
  assert(MS->thePositionInSequence() == 5);
  assert(MS->thePositionInStream()   == 5);
  assert(MS->theSequenceNumber()     == 0);
  assert(strcmp(MS->theFMer().merToString(testmer), "AACTCCGCCCGCACTCTAGC") == 0);

  if (verbose) {
    fprintf(stdout, "MS 6 posInSeq="u64bitFMT" posInStr="u64bitFMT" seqNum="u64bitFMT"\n",
            MS->thePositionInSequence(),
            MS->thePositionInStream(),
            MS->theSequenceNumber());
    if (strcmp(MS->theFMer().merToString(testmer), "AACTCCGCCCGCACTCTAGC"))
      fprintf(stdout, "MS 6 failed '%s' != ''.\n", testmer);
  }

  nm = MS->nextMer();
  assert(nm == false);

  if (verbose && nm)
    fprintf(stdout, "MS 7 failed - still more mers?\n");

  if (verbose)
    fprintf(stdout, "testMSsimple() finished.\n");
}

int
main(int argc, char **argv) {
  mt_s *mtctx = mtInit(451677);

  u32bit nums = 0;
  char **allh = 0L;
  char **alls = 0L;

  if (argc == 0) {
    fprintf(stderr, "usage: %s infile\n", argv[0]);
    exit(1);
  }


  if (0) {
    fprintf(stdout, "seqStream(\"GGGTCAACTCCGCCCGCACTCTAGC\", 25)\n");

    seqStream *SS = new seqStream("GGGTCAACTCCGCCCGCACTCTAGC", 25);

    testSSsimple(SS);
    SS->rewind();
    testSSsimple(SS);

    delete SS;
  }



  if (1) {
    fprintf(stdout, "seqStream(smallfile.fasta)\n");

    FILE *F = fopen("smallfile.fasta", "w");
    fprintf(F, "                  \n");
    fprintf(F, ">sequence1 junk junk\n");
    fprintf(F, "AAAAAAAAAAAA\n");
    fprintf(F, "                  \n");
    fprintf(F, ">sequence2 junk junk\n");
    fprintf(F, "                  \n");
    fprintf(F, "C\nC\nC\nC\nC\nC\nC\nC\nC\nC\nC\nC\n");
    fprintf(F, "                  \n");
    fprintf(F, ">sequence3 junk junk\n");
    fprintf(F, "G G G G G G G G G G G G\n");
    fprintf(F, ">sequence4 junk junk\n");
    fprintf(F, "T  T  T  T  T  T  T  T  T  T  T  T\n");
    fprintf(F, "                  \n");
    fclose(F);

    seqStream *SS = new seqStream("smallfile.fasta");

    testSSsimple(SS);
    SS->rewind();
    testSSsimple(SS);

    delete SS;

    //unlink("smallfile.fasta");
  }



  if (0) {
    fprintf(stdout, "merStream(kMerBuilder(20), seqStream(\"GGGTCAACTCCGCCCGCACTCTAGC\", 25))\n");

    merStream *MS = new merStream(new kMerBuilder(20),
                                  new seqStream("GGGTCAACTCCGCCCGCACTCTAGC", 25),
                                  true, true);

    testMSsimple(MS);
    MS->rewind();
    testMSsimple(MS);
    MS->rewind();
    MS->rewind();
    testMSsimple(MS);

    delete MS;
  }


  exit(0);


  {
    seqFile  *SF = openSeqFile(argv[1]);

    fprintf(stdout, "source '%s' of type '%s' has "u32bitFMT" sequences.\n",
            SF->getSourceName(), SF->getFileTypeName(), SF->getNumberOfSequences());

    fprintf(stdout, "getSequenceLength() vs getSequence(full)\n");
    {
      char  *h = 0L;
      char  *s = 0L;
      u32bit hLen=0, hMax=0;
      u32bit sLen=0, sMax=0;

      for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, h, hLen, hMax, s, sLen, sMax);

        if ((strlen(s) != SF->getSequenceLength(sid)) ||
            (strlen(s) != sLen) ||
            (SF->getSequenceLength(sid) != sLen)) {
          fprintf(stdout, "length differ for sid="u32bitFMT" h='%s' strlen(s)=%d sLen="u32bitFMT" getSequenceLength()="u32bitFMT"\n",
                  sid, h, strlen(s), sLen, SF->getSequenceLength(sid));
        }
      }

      delete [] h;
      delete [] s;
    }


    fprintf(stdout, "getSequenceLength() vs getSequence(part)\n");
    {
      char  *p = new char [128 * 1024 * 1024];

      for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, 0, SF->getSequenceLength(sid), p);

        if (strlen(p) != SF->getSequenceLength(sid)) {
          fprintf(stdout, "length differ for sid="u32bitFMT" strlen(s)=%d getSequenceLength()="u32bitFMT"\n",
                  sid, strlen(p), SF->getSequenceLength(sid));
        }
      }

      delete [] p;
    }


    fprintf(stdout, "loading allh/alls.\n");
    {
      char  *h = 0L;
      char  *s = 0L;
      u32bit hLen=0, hMax=0;
      u32bit sLen=0, sMax=0;

      nums = SF->getNumberOfSequences();
      allh = new char * [SF->getNumberOfSequences()];
      alls = new char * [SF->getNumberOfSequences()];

      for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, h, hLen, hMax, s, sLen, sMax);

        allh[sid] = new char [hLen + 1];
        strcpy(allh[sid], h);

        alls[sid] = new char [sLen + 1];
        strcpy(alls[sid], s);
      }

      delete [] s;
      delete [] h;
    }



    fprintf(stdout, "getSequence(full) vs getSequence(part) & random access\n");
    {
      char  *p = new char [128 * 1024 * 1024];

      for (u32bit tl=0; tl<2 * SF->getNumberOfSequences(); tl++) {
        u32bit sid = mtRandom32(mtctx) % SF->getNumberOfSequences();

        u32bit b   = mtRandom32(mtctx) % SF->getSequenceLength(sid);
        u32bit e   = mtRandom32(mtctx) % SF->getSequenceLength(sid);
        if (b > e) {
          u32bit t=b;
          b = e;
          e = t;
        }

        SF->getSequence(sid, b, e, p);

        if (strlen(p) != e-b) {
          fprintf(stdout, "length differs; strlen(p)=%s requested="u32bitFMT"\n",
                  strlen(p), e-b);
        }

        if (strncmp(alls[sid]+b, p, e-b)) {
          fprintf(stdout, "sequence differs.\n");
        }
      }

      delete [] p;
    }


    delete SF;
  }


  ////////////////////////////////////////


  fprintf(stdout, "seqCache.\n");
  {
    seqCache  *SC = new seqCache(argv[1]);
    delete SC;
  }


  fprintf(stdout, "seqInCore stream.\n");
  {
    seqCache  *SC = new seqCache(argv[1], 0, true);
    seqInCore *IC = SC->getSequenceInCore();

    while (IC) {
      delete IC;
      IC = SC->getSequenceInCore();
    }
  }


  fprintf(stdout, "seqInCore iteration.\n");
  {
    seqCache  *SC = new seqCache(argv[1], 0, true);

    for (u32bit sid=0; sid<SC->getNumberOfSequences(); sid++) {
      seqInCore *IC = SC->getSequenceInCore(sid);
      delete IC;
    }
    delete SC;
  }


  fprintf(stdout, "seqInCore cache random access.\n");
  {
    seqCache  *SC = new seqCache(argv[1], nums/10, true);

    for (u32bit sid=0; sid<SC->getNumberOfSequences(); sid++) {
      u32bit    sid = mtRandom32(mtctx) % SC->getNumberOfSequences();

      //  Cache owns these.
      seqInCore *IC = SC->getSequenceInCore(sid);
    }
    delete SC;
  }


  ////////////////////////////////////////

  fprintf(stdout, "seqStream\n");
  {
    seqStream *S = new seqStream(argv[1]);
    delete S;
  }

  ////////////////////////////////////////

  fprintf(stdout, "merStream\n");
  {
    merStream *M = new merStream(new kMerBuilder(22),
                                 new seqStream(argv[1]),
                                 true, true);
    delete M;
  }

  ////////////////////////////////////////


  for (u32bit sid=0; sid<nums; sid++) {
    delete [] allh[sid];
    delete [] alls[sid];
  }

  delete [] allh;
  delete [] alls;

  return(0);
}

