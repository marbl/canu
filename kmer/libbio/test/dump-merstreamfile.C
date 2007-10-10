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

  seqStream   *STR = new seqStream(argv[1], true);
  seqStore    *STO = new seqStore(argv[1], STR);

#if 0
  fprintf(stdout, "NORMAL\n");
  kMerBuilder *KB1 = new kMerBuilder(8, 0);
#endif
#if 0
  fprintf(stdout, "COMPRESSED\n");
  kMerBuilder *KB1 = new kMerBuilder(8, 1);
#endif
#if 1
  fprintf(stdout, "SPACED\n");  //            01234567890123
  kMerBuilder *KB1 = new kMerBuilder(0, 0, "0011010000");
#endif

  merStream   *MS1 = new merStream(KB1, STO);

  u64bit       s1 = 0;
  u64bit       s2 = 0;


  while (MS1->nextMer()) {
    fprintf(stdout, "'%s' seqNum:"u64bitFMT" posInSeq:"u64bitFMT" posInStream:"u64bitFMT"\n",
            MS1->theFMer().merToString(str),
            MS1->theSequenceNumber(),
            MS1->thePositionInSequence(),
            MS1->thePositionInStream());
    s1 += MS1->theSequenceNumber() + MS1->thePositionInSequence() + MS1->thePositionInStream();
  }

  fprintf(stdout, "REWIND!\n");
  MS1->rewind();

  while (MS1->nextMer()) {
    fprintf(stdout, "'%s' seqNum:"u64bitFMT" posInSeq:"u64bitFMT" posInStream:"u64bitFMT"\n",
            MS1->theFMer().merToString(str),
            MS1->theSequenceNumber(),
            MS1->thePositionInSequence(),
            MS1->thePositionInStream());
    s2 += MS1->theSequenceNumber() + MS1->thePositionInSequence() + MS1->thePositionInStream();
  }

  fprintf(stderr, "s1 = "u64bitFMT"\n", s1);
  fprintf(stderr, "s2 = "u64bitFMT"\n", s2);

  delete STR;
  delete STO;
  delete KB1;
  delete MS1;

  exit(0);
}
