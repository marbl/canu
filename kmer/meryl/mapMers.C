#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libmeryl.H"

int
main(int argc, char **argv) {
  char     *merylFile  = 0L;
  char     *fastaFile  = 0L;
  bool      beVerbose  = false;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merylFile = argv[++arg];
    } else if (strcmp(argv[arg], "-f") == 0) {
      fastaFile = argv[++arg];
    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((merylFile == 0L) || (fastaFile == 0L)) {
    fprintf(stderr, "usage: %s -m meryl -f fasta > output\n", argv[0]);
    exit(1);
  }

  merylStreamReader    *M = new merylStreamReader(merylFile);
  existDB              *E = new existDB(merylFile, M->merSize(), 16);

  FastAWrapper         *F = new FastAWrapper(fastaFile);
  FastASequenceInCore  *S = F->getSequence();

  while (S) {
    merStream  *MS = new merStream(M->merSize(), F->sequence(), F->sequenceLength());

    while (MS->nextMer()) {
      if (E->exists(MS->theFMer())) {
        fprintf(stderr, "%s\tc\t"u32bitFMT"\n", F->defline, MS->thePositionInSequence());
      }
      if (E->exists(MS->theFMer())) {
        fprintf(stderr, "%s\tf\t"u32bitFMT"\n", F->defline, MS->thePositionInSequence());
      }
    }

    delete MS;

    S = F->getSequence();
  }

  delete S;
  delete F;
  delete E;
  delete M;
}
