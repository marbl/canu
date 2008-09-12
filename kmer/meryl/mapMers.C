#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bio++.H"
#include "seqCache.H"
#include "merStream.H"
#include "libmeryl.H"
#include "existDB.H"

#undef FRAGSTATS
#define COVEREDREGIONS

int
main(int argc, char **argv) {
  u32bit    merSize    = 16;
  char     *merylFile  = 0L;
  char     *fastaFile  = 0L;
  bool      beVerbose  = false;
  u32bit    loCount    = 0;
  u32bit    hiCount    = ~u32bitZERO;

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
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((merylFile == 0L) || (fastaFile == 0L)) {
    fprintf(stderr, "usage: %s -m mersize -mers mers -seq fasta > output\n", argv[0]);
    exit(1);
  }

  existDB       *E = new existDB(merylFile, merSize, existDBcounts | existDBcompressCounts | existDBcompressBuckets, loCount, hiCount);
  seqCache      *F = new seqCache(fastaFile);

#ifdef FRAGSTATS
  u32bit         Clen = 0;
  u32bit         Cmax = 4 * 1024 * 1024;
  u32bit        *C    = new u32bit [Cmax];
#endif

  for (u32bit Sid=0; Sid < F->getNumberOfSequences(); Sid++) {
    seqInCore  *S  = F->getSequenceInCore(Sid);
    merStream  *MS = new merStream(new kMerBuilder(merSize),
                                   new seqStream(S->sequence(), S->sequenceLength()),
                                   true, true);

#ifdef COVEREDREGIONS
    u64bit    beg = ~u64bitZERO;
    u64bit    end = ~u64bitZERO;
    u64bit    pos = ~u64bitZERO;

    u64bit    numCovReg = 0;
    u64bit    lenCovReg = 0;
#endif

#ifdef FRAGSTATS
    Clen = 0;
    for (u32bit i=0; i<Cmax; i++)
      C[i] = 0;
#endif

    while (MS->nextMer()) {

#ifdef COVEREDREGIONS
      //  without counts, reports regions with mer coverage.

      //  Orientation tells us nothing, since the mers are probably canonical
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
      }
#endif

#ifdef FRAGSTATS
      //  with counts, report mean, mode, median, min, max for each frag.

      u64bit   cnt = E->count(MS->theFMer()) + E->count(MS->theRMer());

      if (cnt > 0)
        C[Clen++] = cnt;
#endif
    }

#ifdef COVEREDREGIONS
    if (beg != ~u64bitZERO)
      fprintf(stdout, "%s\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\n", S->header(), beg, end+merSize, end+merSize - beg);

    fprintf(stderr, "numCovReg: "u64bitFMT"\n", numCovReg);
    fprintf(stderr, "lenCovReg: "u64bitFMT"\n", lenCovReg);
#endif

#ifdef FRAGSTATS
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
#endif


    delete MS;
    delete S;
  }


  delete F;
  delete E;
}
