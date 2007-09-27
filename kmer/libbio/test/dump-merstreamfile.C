#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bio++.H"

//  Dumps each mer in a merStreamFile.

#define BUILD_SIZE 10
#define TEST_SIZE  10

int
main(int argc, char **argv) {
  char   str[1024];

  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "       dumps mers to stdout.\n");
    exit(1);
  }

  if (0) {
    seqStream   *STR = new seqStream(argv[1], true);
    seqStore    *STO = new seqStore(argv[1], STR);
    kMerBuilder *KB1 = new kMerBuilder(27);
    merStream   *MS1 = new merStream(KB1, STO);

    fprintf(stdout, "NORMAL\n");

    while (MS1->nextMer())
      fprintf(stdout, "'%s' seqNum:"u64bitFMT" posInSeq:"u64bitFMT" posInStream:"u64bitFMT"\n",
              MS1->theFMer().merToString(str),
              MS1->theSequenceNumber(),
              MS1->thePositionInSequence(),
              MS1->thePositionInStream());

    delete STR;
    delete STO;
    delete KB1;
    delete MS1;
  }

  if (0) {
    seqStream   *STR = new seqStream(argv[1], true);
    seqStore    *STO = new seqStore(argv[1], STR);
    kMerBuilder *KB1 = new kMerBuilder(27, 1);
    merStream   *MS1 = new merStream(KB1, STO);

    fprintf(stdout, "COMPRESSED\n");

    while (MS1->nextMer())
      fprintf(stdout, "'%s' seqNum:"u64bitFMT" posInSeq:"u64bitFMT" posInStream:"u64bitFMT"\n",
              MS1->theFMer().merToString(str),
              MS1->theSequenceNumber(),
              MS1->thePositionInSequence(),
              MS1->thePositionInStream());

    delete STR;
    delete STO;
    delete KB1;
    delete MS1;
  }

  if (1) {
    seqStream   *STR = new seqStream(argv[1], true);
    seqStore    *STO = new seqStore(argv[1], STR);
    //                                          01234567890123
    kMerBuilder *KB1 = new kMerBuilder(0, 0, "0011010000");
    merStream   *MS1 = new merStream(KB1, STO);

    fprintf(stdout, "SPACED\n");

    while (MS1->nextMer())
      fprintf(stdout, "'%s' seqNum:"u64bitFMT" posInSeq:"u64bitFMT" posInStream:"u64bitFMT"\n",
              MS1->theFMer().merToString(str),
              MS1->theSequenceNumber(),
              MS1->thePositionInSequence(),
              MS1->thePositionInStream());

    delete STR;
    delete STO;
    delete KB1;
    delete MS1;
  }

  exit(0);
}
