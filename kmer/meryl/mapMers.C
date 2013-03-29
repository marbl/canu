#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bio++.H"
#include "seqCache.H"
#include "merStream.H"
#include "libmeryl.H"
#include "existDB.H"

#define  OP_NONE    0
#define  OP_STATS   1
#define  OP_REGIONS 2
#define  OP_DETAILS 3

int
main(int argc, char **argv) {
  u32bit    merSize    = 16;
  char     *merylFile  = 0L;
  char     *fastaFile  = 0L;
  bool      beVerbose  = false;
  u32bit    loCount    = 0;
  u32bit    hiCount    = ~u32bitZERO;
  u32bit    operation  = OP_NONE;

  //  For OP_STATS

  u32bit         Clen = 0;
  u32bit         Cmax = 4 * 1024 * 1024;
  u32bit        *C    = new u32bit [Cmax];

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merSize = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-mers") == 0) {
      merylFile = argv[++arg];

    } else if (strcmp(argv[arg], "-seq") == 0) {
      fastaFile = argv[++arg];

    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;

    } else if (strcmp(argv[arg], "-lo") == 0) {
      loCount = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-hi") == 0) {
      hiCount = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-stats") == 0) {
      operation = OP_STATS;

    } else if (strcmp(argv[arg], "-regions") == 0) {
      operation = OP_REGIONS;

    } else if (strcmp(argv[arg], "-details") == 0) {
      operation = OP_DETAILS;

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((operation == OP_NONE) || (merylFile == 0L) || (fastaFile == 0L)) {
    fprintf(stderr, "usage: %s [-stats | -regions | -details] -m mersize -mers mers -seq fasta > output\n", argv[0]);
    exit(1);
  }

#if 0
  existDB       *E = NULL;

  if (fileExists("junk.existDB")) {
    fprintf(stderr, "loading from junk.existDB\n");
    E = new existDB("junk.existDB");
    fprintf(stderr, "loaded\n");
  } else {
    exit(1);
    E = new existDB(merylFile, merSize, existDBcounts, loCount, hiCount);
    E->saveState("junk.existDB");
  }
#endif

  existDB  *E = new existDB(merylFile, merSize, existDBcounts, loCount, hiCount);
  seqCache *F = new seqCache(fastaFile);

  fprintf(stderr, "Begin.\n");


  for (u32bit Sid=0; Sid < F->getNumberOfSequences(); Sid++) {
    seqInCore  *S  = F->getSequenceInCore(Sid);
    merStream  *MS = new merStream(new kMerBuilder(merSize),
                                   new seqStream(S->sequence(), S->sequenceLength()),
                                   true, true);

    //  with counts, report mean, mode, median, min, max for each frag.
    if (operation == OP_STATS) {
      Clen = 0;
      for (u32bit i=0; i<Cmax; i++)
        C[i] = 0;

      while (MS->nextMer()) {
        u64bit   cnt = E->count(MS->theFMer()) + E->count(MS->theRMer());

        if (cnt > 0)
          C[Clen++] = cnt;
      }

      u64bit         mean     = u64bitZERO;
      u64bit         min      = ~u64bitZERO;
      u64bit         max      = u64bitZERO;
      u64bit         hist[16]  = { 0 };

      //  Histogram values are powers of two, e.g., <=1, <=2, <=4, <=8, <=16, <=32, <=64, <=128, <=256, <=512, <=1024, <=4096, <=8192, <=328768

      for (u32bit i=0; i<Clen; i++) {
        mean += C[i];

        if ((min > C[i]) && (C[i] > 1))
          min = C[i];
        if (max < C[i])
          max = C[i];

        hist[ logBaseTwo64(C[i]) ]++;
      }

      mean /= Clen;

      fprintf(stdout,
              "%s\t"
              u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"
              u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"
              u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\n",
              S->header(),
              mean, min, max,
              hist[ 0], hist[ 1], hist[ 2], hist[ 3], hist[ 4], hist[ 5], hist[ 6], hist[ 7],
              hist[ 8], hist[ 9], hist[10], hist[11], hist[12], hist[13], hist[14], hist[15]);
    }


    //  without counts, reports regions with mer coverage.
    //  Orientation tells us nothing, since the mers are probably canonical
    if (operation == OP_REGIONS) {
      u64bit    beg = ~u64bitZERO;
      u64bit    end = ~u64bitZERO;
      u64bit    pos = ~u64bitZERO;

      u64bit    numCovReg = 0;
      u64bit    lenCovReg = 0;

      while (MS->nextMer()) {
        if (E->exists(MS->theFMer()) || E->exists(MS->theRMer())) {
          pos = MS->thePositionInSequence();

          if (beg == ~u64bitZERO)
            beg = end = pos;

          if (pos <= end + merSize) {
            end = pos;
          } else {
            fprintf(stdout, "%s\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\n", S->header(), beg, end+merSize, end+merSize - beg);
            numCovReg++;
            lenCovReg += end+merSize - beg;
            beg = end = pos;
          }
        } else {
          fprintf(stdout, "%s\t"u64bitFMT"\tuncovered\n", S->header(), MS->thePositionInSequence());
        }
      }

      if (beg != ~u64bitZERO)
        fprintf(stdout, "%s\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\n", S->header(), beg, end+merSize, end+merSize - beg);

      fprintf(stderr, "numCovReg: "u64bitFMT"\n", numCovReg);
      fprintf(stderr, "lenCovReg: "u64bitFMT"\n", lenCovReg);
    }



    if (operation == OP_DETAILS) {
      char  merString[256];

      while (MS->nextMer()) {
        u64bit beg = MS->thePositionInSequence();
        u64bit end = beg + merSize;
        u64bit fnt = E->count(MS->theFMer());
        u64bit rnt = E->count(MS->theRMer());

        fprintf(stdout, "%s\t%s\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\n",
                S->header(),
                MS->theFMer().merToString(merString),
                beg,
                end,
                fnt,
                rnt,
                fnt + rnt);
      }
    }


    delete MS;
    delete S;
  }


  delete F;
  delete E;
}
