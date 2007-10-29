#include "util++.H"
#include "bio++.H"
#include "existDB.H"

int
main(int argc, char **argv) {
  char  *merFile   = 0L;
  char  *queryFile = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merFile = argv[++arg];
    } else if (strcmp(argv[arg], "-q") == 0) {
      queryFile = argv[++arg];
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
    }
    arg++;
  }

  existDB      *E = new existDB(merFile, 22, existDBnoFlags, 0, ~u32bitZERO);
  seqFile      *Q = openSeqFile(queryFile);
  seqInCore    *S = Q->getSequenceInCore();

  intervalList          IL;
  speedCounter          SC(" %8f frags (%8.5f frags/sec)\r", 1, 1000, true);

  while (S) {
    kMerBuilder    KB(22);
    merStream     *MS = new merStream(&KB, S);

    IL.clear();

    while (MS->nextMer()) {
      if (E->exists(MS->theFMer())) {
        IL.add(MS->thePositionInSequence(), 22);
      }
    }

    IL.merge();

    if (IL.sumOfLengths() > 0) {
      fprintf(stdout, "%5.2f\n",
              100.0 * IL.sumOfLengths() / (double)S->sequenceLength());
    }

    delete MS;
    delete S;

    SC.tick();

    S = Q->getSequenceInCore();
  }

  delete Q;
  delete E;

  return(0);
}

