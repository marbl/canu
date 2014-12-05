#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bio++.H"
#include "seqCache.H"
#include "merStream.H"
#include "libmeryl.H"
#include "existDB.H"

#warning this code might not work due to intervalList changes

int
main(int argc, char **argv) {
  uint32    merSize    = 16;
  char     *merylFile  = 0L;
  char     *fastaFile  = 0L;
  bool      beVerbose  = false;
  uint32    loCount    = 0;
  uint32    hiCount    = ~uint32ZERO;
  uint32    windowsize = 0;
  uint32    skipsize   = 0;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merSize = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-mers") == 0) {
      merylFile = argv[++arg];
    } else if (strcmp(argv[arg], "-seq") == 0) {
      fastaFile = argv[++arg];
    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;
    } else if (strcmp(argv[arg], "-lo") == 0) {
      loCount = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-hi") == 0) {
      hiCount = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-w") == 0) {
      windowsize = strtouint32(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-s") == 0) {
      skipsize = strtouint32(argv[++arg], 0L);
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

  for (uint32 Sid=0; Sid < F->getNumberOfSequences(); Sid++) {
    seqInCore  *S  = F->getSequenceInCore(Sid);
    merStream  *MS = new merStream(new kMerBuilder(merSize),
                                   new seqStream(S->sequence(), S->sequenceLength()),
                                   true, true);

    uint32                         idlen = 0;
    intervalDepthRegions<uint64>  *id    = new intervalDepthRegions<uint64> [S->sequenceLength() * 2 + 2];

    while (MS->nextMer()) {
      int32   cnt = (int32)E->count(MS->theFMer()) + (int32)E->count(MS->theRMer());

      //  Old intervalDepth was to add 'cnt' in the first and subtract 'cnt' in the second.
      //  Then to use the 'ct' field below.
      //  New intervalDepth is the same, but uses the value field.
      //  Count is now the number of intervals that are represented in this block.

      id[idlen].pos     = MS->thePositionInSequence();
      id[idlen].change  = cnt;
      id[idlen].open    = true;
      idlen++;

      id[idlen].pos     = MS->thePositionInSequence() + merSize;
      id[idlen].change  = cnt;
      id[idlen].open    = false;
      idlen++;
    }

    intervalList<uint64>  ID(id, idlen);
    uint32                x = 0;

    uint32                len = S->sequenceLength();

    //  Default case, report un-averaged depth at every single location.
    //
    if ((windowsize == 0) && (skipsize == 0)) {
      for (uint32 i=0; i < ID.numberOfIntervals(); i++) {
        for (; x < ID.lo(i); x++)
          fprintf(stdout, uint32FMTW(7)"\t"uint32FMTW(6)"\n", x, 0);
        for (; x < ID.hi(i); x++)
          fprintf(stdout, uint32FMTW(7)"\t"uint32FMTW(6)"\n", x, ID.value(i));
      }
      for (; x < len; x++)
        fprintf(stdout, uint32FMTW(7)"\t"uint32FMTW(6)"\n", x, 0);

    } else {
      uint32  *depth = new uint32 [len];
      for (x=0; x < len; x++)
        depth[x] = 0;

      for (uint32 i=0; i < ID.numberOfIntervals(); i++)
        for (x=ID.lo(i); x < ID.hi(i); x++)
          depth[x] = ID.count(i);

      uint32   avedepth = 0;

      for (x=0; x < windowsize; x++)
        avedepth += depth[x];

      while (x < len) {
        uint32  avepos = (x - 1) - (windowsize - 1) / 2;
        if ((avepos % skipsize) == 0)
          fprintf(stdout, uint32FMT"\t%.4f\n",
                  avepos,
                  (double)avedepth / (double)windowsize);

        avedepth = avedepth + depth[x] - depth[x-windowsize];

        x++;
      }

      delete [] depth;
    }

    delete [] id;

    delete MS;
    delete S;
  }


  delete F;
  delete E;
}
