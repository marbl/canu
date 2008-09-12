#include "bio++.H"
#include "seqStream.H"

int
main(int argc, char **argv) {
  seqStream    A;
  seqStream   *B = 0L;
  seqStream   *C = 0L;

  seqCache    *F = new seqCache(argv[1]);

  B = new seqStream(argv[1], true);
  C = new seqStream(F, true);

  A.setFile(argv[1]);
  A.parse("0,2,1");
  A.finish();

  A.setSeparator('.', 10);

  fprintf(stderr, "A: numSeqs "u32bitFMT" len "u64bitFMT"\n", A.numberOfSequences(), A.lengthOfSequences());
  fprintf(stderr, "B: numSeqs "u32bitFMT" len "u64bitFMT"\n", B->numberOfSequences(), B->lengthOfSequences());
  fprintf(stderr, "C: numSeqs "u32bitFMT" len "u64bitFMT"\n", C->numberOfSequences(), C->lengthOfSequences());

  for (u32bit i=0; i<A.numberOfSequences(); i++) {
    printf("seq:"u32bitFMT" len:"u64bitFMT" start:"u64bitFMT" iid:"u64bitFMT"\n",
           i, A.lengthOf(i), A.startOf(i), A.IIDOf(i));
  }

  char   ch = 0;
  u32bit xx = 0;
  while (ch = A.get()) {
    u64bit sp = A.seqPos();
    u64bit si = A.seqIID();
    u64bit st = A.strPos();
    printf("%c -- seqPos:"u64bitFMT" seqIID:"u64bitFMT" strPos:"u64bitFMT" lookup:"u32bitFMT"\n",
           ch, sp, si, st,
           A.sequenceNumberOfPosition(xx));
    xx++;
  }

  A.rewind();

  while (!A.eof()) {
    printf("%c", A.get());
  }
  printf("\n");

  A.rewind();

  //  24,3
  kMerBuilder  KB(8, 3);
  merStream   *MS = new merStream(&KB, &A);

  while (MS->nextMer()) {
    char  mstr[1024];

    u64bit sp = MS->thePositionInSequence();
    u64bit si = MS->theSequenceNumber();
    u64bit st = MS->thePositionInStream();

    printf("%s -- seqPos:"u64bitFMT" seqIID:"u64bitFMT" strPos:"u64bitFMT"\n",
           MS->theFMer().merToString(mstr),
           sp, si, st);

  }

  delete C;
  delete B;

  return(0);
}
