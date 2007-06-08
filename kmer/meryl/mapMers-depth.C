#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bio++.H"
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

  existDB       *E = new existDB(merylFile, merSize, 22, loCount, hiCount, existDBcounts | existDBcompressCounts | existDBcompressBuckets);
  seqFile       *F = openSeqFile(fastaFile);
  seqInCore     *S = F->getSequenceInCore();

  char           merstringf[256] = {0};
  char           merstringr[256] = {0};

  while (S) {
    merStream             *MS = new merStream(merSize, S);

    u32bit                 idlen = 0;
    intervalDepthRegions  *id    = new intervalDepthRegions [S->sequenceLength() * 2 + 2];

    while (MS->nextMer()) {
      s32bit   cnt = (s32bit)E->count(MS->theFMer()) + (s32bit)E->count(MS->theRMer());

#if 0
      if (idlen < 150)
        fprintf(stdout, "add %d %d %d %d %d %s %s\n",
                (int)idlen,
                (int)MS->thePositionInSequence(),
                (int)cnt,
                E->count(MS->theFMer()),
                E->count(MS->theRMer()),
                MS->theFMer().merToString(merstringf),
                MS->theRMer().merToString(merstringr));
#endif

      id[idlen].pos = MS->thePositionInSequence();
      id[idlen].cha = cnt;
      idlen++;

      id[idlen].pos = MS->thePositionInSequence() + merSize;
      id[idlen].cha = -cnt;
      idlen++;
    }

    intervalDepth ID(id, idlen);
    for (u32bit i=0; i<ID.numberOfIntervals(); i++)
      for (u32bit x=ID.lo(i); x<ID.hi(i); x++)
        //fprintf(stdout, u32bitFMTW(7)"\t"u64bitFMTW(7)"\t"u64bitFMTW(7)"\t"u32bitFMTW(6)"\n", i, ID.lo(i), ID.hi(i), ID.de(i));
        fprintf(stdout, u32bitFMTW(7)"\t"u32bitFMTW(6)"\n", x, ID.de(i));

    delete [] id;

    delete MS;
    delete S;

    S = F->getSequenceInCore();
  }

  delete S;
  delete F;
  delete E;
}
