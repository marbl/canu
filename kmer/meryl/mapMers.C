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
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((merylFile == 0L) || (fastaFile == 0L)) {
    fprintf(stderr, "usage: %s -m mersize -mers mers -seq fasta > output\n", argv[0]);
    exit(1);
  }

  char                  merstring[1024];

  existDB              *E = new existDB(merylFile, merSize, 16);
  FastABase            *F = new FastAFile(fastaFile);
  FastASequenceInCore  *S = F->getSequence();

  while (S) {
    merStream  *MS = new merStream(merSize, S->sequence(), S->sequenceLength());

    u64bit    beg = ~u64bitZERO;
    u64bit    end = ~u64bitZERO;
    u64bit    pos = ~u64bitZERO;

    while (MS->nextMer()) {

      //  Orientation tells us nothing, since the mers are probably canonical

      if (E->exists(MS->theFMer()) || E->exists(MS->theRMer())) {
        pos = MS->thePositionInSequence();

        if (beg == ~u64bitZERO)
          beg = end = pos;

        if (pos <= end + merSize) {
          end = pos;
        } else {
          fprintf(stdout, "%s\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\n", S->header(), beg, end+22, end+22 - beg);
          beg = end = pos;
        }
      }
    }

    if (beg != ~u64bitZERO)
      fprintf(stdout, "%s\t"u64bitFMT"\t"u64bitFMT"\t"u64bitFMT"\n", S->header(), beg, end+22, end+22 - beg);

    delete MS;
    delete S;

    S = F->getSequence();
  }

  delete S;
  delete F;
  delete E;
}
