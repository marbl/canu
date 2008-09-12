#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bio++.H"
#include "seqCache.H"
#include "merStream.H"
#include "libmeryl.H"
#include "existDB.H"


int
main(int argc, char **argv) {
  u32bit    merSize    = 16;
  char     *merylFile  = 0L;
  char     *fastaFile  = 0L;
  bool      beVerbose  = false;
  u32bit    loCount    = 0;
  u32bit    hiCount    = ~u32bitZERO;
  u32bit    windowsize = 0;
  u32bit    skipsize   = 0;

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
    } else if (strcmp(argv[arg], "-w") == 0) {
      windowsize = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-s") == 0) {
      skipsize = strtou32bit(argv[++arg], 0L);
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

  for (u32bit Sid=0; Sid < F->getNumberOfSequences(); Sid++) {
    seqInCore  *S  = F->getSequenceInCore(Sid);
    merStream  *MS = new merStream(new kMerBuilder(merSize),
                                   new seqStream(S->sequence(), S->sequenceLength()),
                                   true, true);

    u32bit                 idlen = 0;
    intervalDepthRegions  *id    = new intervalDepthRegions [S->sequenceLength() * 2 + 2];

    while (MS->nextMer()) {
      s32bit   cnt = (s32bit)E->count(MS->theFMer()) + (s32bit)E->count(MS->theRMer());

      id[idlen].pos = MS->thePositionInSequence();
      id[idlen].cha = cnt;
      idlen++;

      id[idlen].pos = MS->thePositionInSequence() + merSize;
      id[idlen].cha = -cnt;
      idlen++;
    }

    intervalDepth ID(id, idlen);
    u32bit        x = 0;

    u32bit        len = S->sequenceLength();

    //  Default case, report un-averaged depth at every single location.
    //
    if ((windowsize == 0) && (skipsize == 0)) {
      for (u32bit i=0; i < ID.numberOfIntervals(); i++) {
        for (; x < ID.lo(i); x++)
          fprintf(stdout, u32bitFMTW(7)"\t"u32bitFMTW(6)"\n", x, 0);
        for (; x < ID.hi(i); x++)
          fprintf(stdout, u32bitFMTW(7)"\t"u32bitFMTW(6)"\n", x, ID.de(i));
      }
      for (; x < len; x++)
        fprintf(stdout, u32bitFMTW(7)"\t"u32bitFMTW(6)"\n", x, 0);

    } else {
      u32bit  *depth = new u32bit [len];
      for (x=0; x < len; x++)
        depth[x] = 0;

      for (u32bit i=0; i < ID.numberOfIntervals(); i++)
        for (x=ID.lo(i); x < ID.hi(i); x++)
          depth[x] = ID.de(i);

      u32bit   avedepth = 0;

      for (x=0; x < windowsize; x++)
        avedepth += depth[x];

      while (x < len) {
        u32bit  avepos = (x - 1) - (windowsize - 1) / 2;
        if ((avepos % skipsize) == 0)
          fprintf(stdout, u32bitFMT"\t%.4f\n",
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
